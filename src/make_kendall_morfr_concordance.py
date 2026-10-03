import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from scipy.stats import kendalltau, hypergeom
from statsmodels.stats.multitest import multipletests


TARGETS = {
	'Age': 'Age',
	'SjD Diagnosis': 'SjD_Diagnosis',
	'Focus Score': 'Focus_Score',
	'SSA': 'SSA',
	'SSB': 'SSB',
	'ANA': 'ANA',
	'RF': 'RF',
}
TOP_K = [25, 50, 100, 250, 500, 1000]
TARGET_COLORS = {
	'Age': '#397a78',
	'SjD Diagnosis': '#c36f4d',
	'Focus Score': '#5978a8',
	'SSA': '#9b6a40',
	'SSB': '#6a8f53',
	'ANA': '#8a617e',
	'RF': '#687078',
}


def save_figure(figure, stem):
	figure.savefig(stem.with_suffix('.pdf'), bbox_inches='tight')
	figure.savefig(stem.with_suffix('.png'), dpi=300, bbox_inches='tight')
	plt.close(figure)


def calculate_concordance(var1_path, morfr_path, counts_path):
	kendall_scores = pd.read_csv(var1_path, index_col=0)
	kendall_scores.index = kendall_scores.index.astype(str)
	counts = pd.read_csv(counts_path, index_col=0)
	gene_universe = set(counts.index.astype(str))
	morfr = pd.read_csv(morfr_path)
	morfr['Gene'] = morfr['Gene'].astype(str)
	overlap_rows = []
	tau_rows = []
	for target, score_column in TARGETS.items():
		scores = pd.to_numeric(kendall_scores[score_column], errors='coerce').dropna()
		common_genes = set(scores.index).intersection(gene_universe)
		scores = scores.loc[sorted(common_genes)].sort_values(ascending=False, kind='mergesort')
		tree_rankings = morfr[morfr['Target'] == target].sort_values('MORFR_Rank')
		for k in TOP_K:
			kendall_top = set(scores.head(k).index)
			tree_top = tree_rankings.head(k)['Gene'].tolist()
			tree_top_common = set(tree_top).intersection(common_genes)
			observed = len(kendall_top.intersection(tree_top_common))
			expected = k * len(tree_top_common) / len(common_genes)
			pvalue = float(hypergeom.sf(observed - 1, len(common_genes), len(tree_top_common), k))
			overlap_rows.append({
				'Metric': 'Top-k overlap',
				'Target': target,
				'K': k,
				'Kendall_Universe_N': len(common_genes),
				'MORFR_TopK_In_Universe': len(tree_top_common),
				'Observed_Overlap': observed,
				'Expected_Overlap_Random': expected,
				'Fold_Enrichment': observed / expected if expected else np.nan,
				'P_Hypergeometric': pvalue,
			})
		matched_tree = tree_rankings[tree_rankings['Gene'].isin(common_genes)].copy()
		matched_tree['Kendall_Score'] = matched_tree['Gene'].map(scores)
		tau, pvalue = kendalltau(-matched_tree['MORFR_Rank'], matched_tree['Kendall_Score'])
		tau_rows.append({
			'Metric': 'Kendall tau within MORFR top 1,000',
			'Target': target,
			'K': 1000,
			'Matched_Genes_N': len(matched_tree),
			'Kendall_Tau': float(tau),
			'P_Kendall_Tau': float(pvalue),
		})

	overlap = pd.DataFrame(overlap_rows)
	overlap['Q_BH_42_TopK_Tests'] = multipletests(overlap['P_Hypergeometric'], method='fdr_bh')[1]
	tau_results = pd.DataFrame(tau_rows)
	tau_results['Q_BH_7_Tau_Tests'] = multipletests(tau_results['P_Kendall_Tau'], method='fdr_bh')[1]
	return overlap, tau_results


def plot_concordance(overlap, tau_results, output_stem):
	target_order = list(TARGETS)
	k_order = TOP_K
	figure, (heat_axis, tau_axis) = plt.subplots(
		1, 2, figsize=(12.4, 6.0),
		gridspec_kw={'width_ratios': [1.35, 1.0]},
	)
	lift_matrix = np.zeros((len(target_order), len(k_order)))
	for row_index, target in enumerate(target_order):
		for column_index, k in enumerate(k_order):
			row = overlap[(overlap['Target'] == target) & (overlap['K'] == k)].iloc[0]
			lift_matrix[row_index, column_index] = np.log2(max(row['Fold_Enrichment'], 1e-12))
	image = heat_axis.imshow(lift_matrix, cmap='YlGnBu', norm=Normalize(0, max(1, float(lift_matrix.max()))), aspect='auto')
	heat_axis.set_xticks(np.arange(len(k_order)), [str(k) for k in k_order])
	heat_axis.set_yticks(np.arange(len(target_order)), target_order)
	heat_axis.set_xlabel('Top-k genes in each ranking')
	heat_axis.set_title('Top-k overlap above random expectation', loc='left', fontweight='bold')
	for row_index, target in enumerate(target_order):
		for column_index, k in enumerate(k_order):
			row = overlap[(overlap['Target'] == target) & (overlap['K'] == k)].iloc[0]
			expected = row['Expected_Overlap_Random']
			expected_label = f'{expected:.2g}' if expected < 1 else (f'{expected:.1f}' if expected < 10 else f'{expected:.0f}')
			label = f"{int(row['Observed_Overlap'])}/{expected_label}"
			color = 'white' if lift_matrix[row_index, column_index] > lift_matrix.max() * 0.55 else '#222222'
			heat_axis.text(column_index, row_index, label, ha='center', va='center', fontsize=6.2, color=color)
	colorbar = figure.colorbar(image, ax=heat_axis, fraction=0.046, pad=0.04)
	colorbar.set_label('log2(fold enrichment)')
	heat_axis.spines[['top', 'right', 'left', 'bottom']].set_visible(False)

	tau_results = tau_results.set_index('Target').loc[target_order].reset_index()
	y_positions = np.arange(len(target_order))
	for y, (_, row) in zip(y_positions, tau_results.iterrows()):
		tau_axis.scatter(row['Kendall_Tau'], y, s=46, color=TARGET_COLORS[row['Target']], zorder=3)
		tau_axis.text(row['Kendall_Tau'] + 0.012, y, f"tau={row['Kendall_Tau']:.2f}; q={row['Q_BH_7_Tau_Tests']:.2g}", va='center', fontsize=7)
	tau_axis.axvline(0, color='#777777', linewidth=0.7)
	tau_axis.set_xlim(-0.02, max(0.40, float(tau_results['Kendall_Tau'].max()) + 0.14))
	tau_axis.set_yticks(y_positions, target_order)
	tau_axis.invert_yaxis()
	tau_axis.set_xlabel('Tie-aware Kendall tau-b (higher score vs better MORFR rank)')
	tau_axis.set_title('Within the saved MORFR top 1,000', loc='left', fontweight='bold')
	tau_axis.grid(axis='x', color='#eeeeee', linewidth=0.6)
	tau_axis.spines[['top', 'right', 'left']].set_visible(False)

	figure.suptitle('Kendall rank scores and MORFR split rankings', x=0.04, ha='left', y=0.99, fontsize=11, fontweight='bold')
	figure.subplots_adjust(left=0.12, right=0.98, top=0.90, bottom=0.10, wspace=0.42)
	save_figure(figure, output_stem)


def main():
	parser = argparse.ArgumentParser(description='Compare all-sample Kendall-rho-squared ranks to fold-safe MORFR split rankings.')
	parser.add_argument('--var1', type=Path, default=Path('var1.csv'))
	parser.add_argument('--morfr', type=Path, default=Path('coding_exon_snp_scan_20261002/morfr_rankings_scored_20261002/morfr_top1000_by_target.csv'))
	parser.add_argument('--counts', type=Path, default=Path('all_gene.csv'))
	parser.add_argument('--table', type=Path, default=Path('coding_exon_sequence_allele_scan_20261002/kendall_morfr_concordance.csv'))
	parser.add_argument('--figure-dir', type=Path, default=Path('figures/publication_20261002'))
	parser.add_argument('--overwrite', action='store_true')
	args = parser.parse_args()
	for path in [args.var1, args.morfr, args.counts]:
		if not path.is_file():
			parser.error(f'Input file does not exist: {path}')
	output_stem = args.figure_dir / 'kendall_morfr_rank_concordance'
	outputs = [output_stem.with_suffix('.pdf'), output_stem.with_suffix('.png'), args.table]
	if not args.overwrite and any(path.exists() for path in outputs):
		parser.error('Concordance output exists; pass --overwrite to replace it')
	args.figure_dir.mkdir(parents=True, exist_ok=True)
	matplotlib.rcParams.update({
		'font.family': 'DejaVu Sans', 'font.size': 8,
		'axes.labelsize': 8, 'axes.titlesize': 9,
		'pdf.fonttype': 42, 'ps.fonttype': 42,
		'axes.linewidth': 0.8, 'savefig.facecolor': 'white',
	})
	overlap, tau_results = calculate_concordance(args.var1, args.morfr, args.counts)
	results = pd.concat([overlap, tau_results], ignore_index=True, sort=False)
	args.table.parent.mkdir(parents=True, exist_ok=True)
	results.to_csv(args.table, index=False)
	plot_concordance(overlap, tau_results, output_stem)
	print(f'Wrote rank-concordance table: {args.table}')
	print(f'Wrote {output_stem.with_suffix(".pdf")} and {output_stem.with_suffix(".png")}')
	print('Within-top-1,000 Kendall tau range:', f"{tau_results['Kendall_Tau'].min():.3f} to {tau_results['Kendall_Tau'].max():.3f}")
	print('BH-adjusted top-k overlap q<0.05 tests:', int((overlap['Q_BH_42_TopK_Tests'] < 0.05).sum()), 'of', len(overlap))


if __name__ == '__main__':
	main()
