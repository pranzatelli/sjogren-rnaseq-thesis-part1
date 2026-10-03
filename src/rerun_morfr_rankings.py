import argparse
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelEncoder


TARGETS = [
	('Age', 'Age', True),
	('SjD Diagnosis', 'SjD_Diagnosis', False),
	('Focus Score', 'Focus_Score', True),
	('SSA', 'SSA', False),
	('SSB', 'SSB', False),
	('ANA', 'ANA', False),
	('RF', 'RF', False),
]


def load_morfr(path):
	spec = importlib.util.spec_from_file_location('swarMORFR_module', path)
	module = importlib.util.module_from_spec(spec)
	assert spec.loader is not None
	spec.loader.exec_module(module)
	return module


def encode_target(target, regression, encoder):
	if regression:
		return pd.to_numeric(target, errors='coerce')
	return pd.Series(encoder.fit_transform(target), index=target.index, name=target.name)


def main():
	parser = argparse.ArgumentParser(
		description='Rerun fold-local MORFR decision-tree feature ranking and intersect with q-tier variant genes.'
	)
	parser.add_argument('--script', type=Path, default=Path('swarMORFR.py'))
	parser.add_argument('--metadata', type=Path, default=Path('meta.csv'))
	parser.add_argument('--gene-counts', type=Path, default=Path('all_gene.csv'))
	parser.add_argument('--candidate-evidence', type=Path, default=Path('coding_exon_snp_scan_20261002/exploratory_evidence_q25.csv'))
	parser.add_argument('--repeats', type=int, default=5)
	parser.add_argument('--trees', type=int, default=100)
	parser.add_argument('--features', type=int, default=1000)
	parser.add_argument('--targets', nargs='+', choices=[target[0] for target in TARGETS], default=[target[0] for target in TARGETS])
	parser.add_argument('--seed', type=int, default=20261002)
	parser.add_argument('--outdir', type=Path, default=Path('coding_exon_snp_scan_20261002'))
	args = parser.parse_args()
	for path in (args.script, args.metadata, args.gene_counts, args.candidate_evidence):
		if not path.is_file():
			parser.error(f'Input file does not exist: {path}')
	if args.repeats < 1 or args.trees < 1 or args.features < 1:
		parser.error('--repeats, --trees and --features must be positive')

	full_path = args.outdir / 'morfr_top1000_by_target.csv'
	candidate_path = args.outdir / 'candidate_morfr_evidence.csv'
	performance_path = args.outdir / 'morfr_target_performance.csv'
	if full_path.exists() or candidate_path.exists() or performance_path.exists():
		parser.error('MORFR outputs already exist; choose a new --outdir to avoid overwrite')

	morfr = load_morfr(args.script)
	metadata = pd.read_csv(args.metadata, index_col=0)
	candidate_evidence = pd.read_csv(args.candidate_evidence)
	candidate_genes = sorted(set(candidate_evidence['Gene'].dropna().astype(str)) - {''})
	ranking_frames = []
	performance_rows = []

	selected_targets = [target for target in TARGETS if target[0] in args.targets]
	for target_index, (target_name, cli_name, regression) in enumerate(selected_targets):
		features, target, diagnosis = morfr.make_data(target_name, enst=False, lim=False)
		target = encode_target(target, regression, LabelEncoder())
		valid = target.notna()
		features = features.loc[valid]
		target = target.loc[valid]
		diagnosis = diagnosis.loc[valid]
		feature_names = features.columns
		importance_sum = np.zeros(len(feature_names), dtype=float)
		occurrence_sum = np.zeros(len(feature_names), dtype=float)
		holdout_scores = []
		for repeat in range(args.repeats):
			np.random.seed(args.seed + target_index * args.repeats + repeat)
			importance, occurrences, holdout_score = morfr.build_model(
				features,
				target,
				diagnosis,
				reg=regression,
				ntrees=args.trees,
				enriched=True,
				balance=True,
				features_per_tree=args.features,
				return_test_score=True,
			)
			importance_sum += importance
			occurrence_sum += occurrences
			holdout_scores.append(holdout_score)
		mean_importance = np.divide(
			importance_sum,
			occurrence_sum,
			out=np.zeros_like(importance_sum),
			where=occurrence_sum > 0,
		)
		descriptive_rho2 = features.corrwith(target, method='kendall').pow(2).reindex(feature_names)
		ranking = pd.DataFrame({
			'Gene': feature_names,
			'Target': target_name,
			'MORFR_Importance': mean_importance,
			'MORFR_Occurrences': occurrence_sum,
			'Kendall_rho2_all_samples_descriptive': descriptive_rho2.to_numpy(),
			'N_samples': len(target),
			'Repeats': args.repeats,
			'Trees_per_repeat': args.trees,
			'Features_sampled_per_tree': min(args.features, len(feature_names)),
		})
		ranking = ranking.sort_values('MORFR_Importance', ascending=False)
		ranking['MORFR_Rank'] = np.arange(1, len(ranking) + 1)
		ranking_frames.append(ranking.head(1000))
		performance_rows.append({
			'Target': target_name,
			'N_samples': len(target),
			'Repeats': args.repeats,
			'Trees_per_repeat': args.trees,
			'Outer_holdout_corr2_mean': float(np.mean(holdout_scores)),
			'Outer_holdout_corr2_sd': float(np.std(holdout_scores, ddof=1)) if len(holdout_scores) > 1 else 0.0,
			'Outer_holdout_corr2_by_repeat': ';'.join(f'{score:.4f}' for score in holdout_scores),
		})
		print(f'{target_name}: ranked {len(ranking)} genes over {args.repeats} x {args.trees} trees', flush=True)

	all_rankings = pd.concat(ranking_frames, ignore_index=True)
	performance = pd.DataFrame(performance_rows)
	args.outdir.mkdir(parents=True, exist_ok=True)
	all_rankings.to_csv(full_path, index=False)
	performance.to_csv(performance_path, index=False)
	joined = candidate_evidence.merge(
		all_rankings,
		on='Gene',
		how='left',
		suffixes=('', '_MORFR'),
	)
	joined = joined.merge(performance, on='Target', how='left')
	joined.to_csv(candidate_path, index=False)
	print(f'Wrote top gene rankings: {full_path}')
	print(f'Wrote held-out performance summary: {performance_path}')
	print(f'Wrote variant-gene/MORFR cross-evidence: {candidate_path}')
	print('MORFR importance is exploratory and based on a single outer holdout per repeat; Kendall rho^2 is descriptive and uses all samples.')


if __name__ == '__main__':
	main()