import argparse
from pathlib import Path

import pandas as pd


PHENOTYPES = ['SjD Diagnosis', 'SSA', 'SSB', 'ANA', 'RF']
VARIANT_COLUMNS = [
	'CHROM', 'POS', 'REF', 'ALT', 'VARIANT_TYPE', 'INDEL_LENGTH',
	'GENES', 'N_CARRIERS', 'N_HOM_REF', 'CALL_RATE',
]


def tier_for_q(value):
	if value < 0.05:
		return 'Strong (q<0.05)'
	if value < 0.10:
		return 'Suggestive (0.05<=q<0.10)'
	return 'Exploratory (0.10<=q<0.25)'


def main():
	parser = argparse.ArgumentParser(
		description='Join pooled variant q tiers to existing held-out MORFR rankings.'
	)
	parser.add_argument(
		'--scan', type=Path,
		default=Path('coding_exon_sequence_allele_scan_20261002/all_tests.csv.gz'),
	)
	parser.add_argument(
		'--morfr-dir', type=Path,
		default=Path('coding_exon_snp_scan_20261002/morfr_rankings_scored_20261002'),
	)
	parser.add_argument(
		'--out', type=Path,
		default=Path('coding_exon_sequence_allele_scan_20261002/morfr_tier_evidence.csv'),
	)
	args = parser.parse_args()
	inputs = [
		args.scan,
		args.morfr_dir / 'morfr_top1000_by_target.csv',
		args.morfr_dir / 'morfr_target_performance.csv',
	]
	for path in inputs:
		if not path.is_file():
			parser.error(f'Input file does not exist: {path}')
	if args.out.exists():
		parser.error(f'Output already exists; choose a new path: {args.out}')

	scan = pd.read_csv(args.scan)
	q_columns = [f'{phenotype}_Q_BH_ALL_TESTS' for phenotype in PHENOTYPES]
	missing = [column for column in VARIANT_COLUMNS + q_columns if column not in scan]
	if missing:
		parser.error('Scan is missing columns: ' + ', '.join(missing))
	associations = scan.melt(
		id_vars=VARIANT_COLUMNS,
		value_vars=q_columns,
		var_name='Phenotype_Q',
		value_name='Variant_Q_BH_All_Tests',
	)
	associations = associations[associations['Variant_Q_BH_All_Tests'] < 0.25].copy()
	associations['Target'] = associations['Phenotype_Q'].str.replace(
		'_Q_BH_ALL_TESTS', '', regex=False
	)
	associations['Variant_Tier'] = associations['Variant_Q_BH_All_Tests'].map(tier_for_q)
	associations['Gene'] = associations['GENES'].str.split(';')
	associations = associations.explode('Gene')

	rankings = pd.read_csv(inputs[1])
	performance = pd.read_csv(inputs[2])
	associations = associations.merge(
		rankings[['Gene', 'Target', 'MORFR_Rank', 'MORFR_Importance', 'MORFR_Occurrences']],
		on=['Gene', 'Target'],
		how='left',
	)
	gene_top1000 = (
		rankings.sort_values('MORFR_Rank')
		.groupby('Gene', as_index=False)
		.agg(
			MORFR_Top1000_Target_Count=('Target', 'nunique'),
			MORFR_Best_Top1000_Target=('Target', 'first'),
			MORFR_Best_Top1000_Rank=('MORFR_Rank', 'first'),
		)
	)
	associations = associations.merge(gene_top1000, on='Gene', how='left')
	associations = associations.merge(
		performance[['Target', 'Outer_holdout_corr2_mean', 'Outer_holdout_corr2_sd']],
		on='Target',
		how='left',
	)
	associations = associations.sort_values(
		['Variant_Q_BH_All_Tests', 'CHROM', 'POS', 'Gene']
	).reset_index(drop=True)
	args.out.parent.mkdir(parents=True, exist_ok=True)
	associations.to_csv(args.out, index=False)

	unique_tests = associations.drop_duplicates(
		['CHROM', 'POS', 'REF', 'ALT', 'Target', 'Variant_Q_BH_All_Tests']
	)
	print(f'Wrote {len(associations)} variant-gene rows for {len(unique_tests)} variant-phenotype tests to {args.out}')
	for tier, group in unique_tests.groupby('Variant_Tier', sort=False):
		keys = group[['CHROM', 'POS', 'REF', 'ALT', 'Target']].drop_duplicates()
		gene_rows = associations.merge(keys, on=['CHROM', 'POS', 'REF', 'ALT', 'Target'])
		matched = int(gene_rows['MORFR_Rank'].notna().sum())
		print(f'{tier}: {len(keys)} tests; {matched}/{len(gene_rows)} variant-gene rows in matching-target MORFR top 1,000')


if __name__ == '__main__':
	main()