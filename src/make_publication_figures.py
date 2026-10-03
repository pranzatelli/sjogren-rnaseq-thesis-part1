import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


PHENOTYPES = ['SjD Diagnosis', 'SSA', 'SSB', 'ANA', 'RF']
DE_FILES = {
	'SjD diagnosis': 'gene.disease.csv',
	'SSA': 'gene.ssa.csv',
	'SSB': 'gene.ssb.csv',
	'ANA': 'gene.ana.csv',
	'RF': 'gene.rf.csv',
}
TIER_COLORS = {
	'Strong (q<0.05)': '#b4474b',
	'Suggestive (0.05<=q<0.10)': '#d28b19',
	'Exploratory (0.10<=q<0.25)': '#3e7884',
}
VARIANT_MARKERS = {'SNV': 'o', 'MNV': 'D', 'insertion': '^', 'deletion': 'v'}


def tier_for_q(value):
	if value < 0.05:
		return 'Strong (q<0.05)'
	if value < 0.10:
		return 'Suggestive (0.05<=q<0.10)'
	return 'Exploratory (0.10<=q<0.25)'


def save_figure(figure, stem):
	figure.savefig(stem.with_suffix('.pdf'), bbox_inches='tight')
	figure.savefig(stem.with_suffix('.png'), dpi=300, bbox_inches='tight')
	plt.close(figure)


def load_variant_associations(path):
	scan = pd.read_csv(path)
	q_columns = [f'{phenotype}_Q_BH_ALL_TESTS' for phenotype in PHENOTYPES]
	missing = [column for column in q_columns if column not in scan]
	if missing:
		raise ValueError('Variant scan lacks pooled q columns: ' + ', '.join(missing))
	id_columns = [
		'CHROM', 'POS', 'REF', 'ALT', 'VARIANT_TYPE', 'INDEL_LENGTH',
		'GENES', 'N_CARRIERS', 'N_HOM_REF', 'CALL_RATE',
	]
	associations = scan.melt(
		id_vars=id_columns,
		value_vars=q_columns,
		var_name='Phenotype_Q',
		value_name='q',
	)
	associations = associations[associations['q'] < 0.25].copy()
	associations['Phenotype'] = associations['Phenotype_Q'].str.replace(
		'_Q_BH_ALL_TESTS', '', regex=False
	)
	associations['Tier'] = associations['q'].map(tier_for_q)
	associations['Gene_List'] = associations['GENES'].fillna('').str.split(';')
	associations['Display_Gene'] = associations['Gene_List'].map(
		lambda genes: genes[0] + (f' +{len(genes) - 1}' if len(genes) > 1 else '')
	)
	associations['Label'] = associations.apply(
		lambda row: f"{row['Display_Gene']} / {row['Phenotype']} / {row['CHROM']}:{row['POS']}",
		axis=1,
	)
	return associations


def plot_variant_tiers(associations, path):
	associations = associations.sort_values('q', ascending=True).reset_index(drop=True)
	associations['neg_log10_q'] = -np.log10(associations['q'].clip(lower=1e-300))
	associations['row'] = np.arange(len(associations))
	figure, axis = plt.subplots(figsize=(10.5, max(5.8, len(associations) * 0.34)))
	for _, row in associations.iterrows():
		axis.scatter(
			row['neg_log10_q'], row['row'],
			color=TIER_COLORS[row['Tier']],
			marker=VARIANT_MARKERS[row['VARIANT_TYPE']],
			s=42, linewidth=0.7, edgecolor='white', zorder=3,
		)
		axis.text(row['neg_log10_q'] + 0.055, row['row'], f"q={row['q']:.3g}", va='center', fontsize=7.5)

	for threshold in (0.25, 0.10, 0.05):
		axis.axvline(-np.log10(threshold), color='#777777', linestyle='--', linewidth=0.8, zorder=0)

	axis.set_yticks(associations['row'])
	axis.set_yticklabels(associations['Label'], fontsize=8)
	axis.invert_yaxis()
	axis.set_xlabel('-log10 pooled BH-adjusted CMH q-value')
	axis.set_ylabel('Variant-gene / phenotype')
	axis.set_title('Exon-overlapping sequence-allele associations', loc='left', fontweight='bold', pad=16)
	axis.grid(axis='y', color='#e6e6e6', linewidth=0.5)
	axis.grid(axis='x', color='#eeeeee', linewidth=0.5)
	axis.spines[['top', 'right']].set_visible(False)

	legend = [
		Line2D([0], [0], marker='o', color='none', markerfacecolor=TIER_COLORS[tier],
		       markeredgecolor='white', markersize=7, label=label)
		for tier, label in [
			('Strong (q<0.05)', '<0.05'),
			('Suggestive (0.05<=q<0.10)', '<0.1'),
		]
	]
	legend.extend(
		Line2D([0], [0], marker=marker, color='#555555', linestyle='none',
		       markersize=6, label=variant_type)
		for variant_type, marker in VARIANT_MARKERS.items()
		if variant_type in set(associations['VARIANT_TYPE'])
	)
	axis.legend(handles=legend, frameon=False, ncol=4, loc='upper center',
	            bbox_to_anchor=(0.52, -0.08), fontsize=7.5)
	figure.subplots_adjust(left=0.34, right=0.94, top=0.91, bottom=0.14)
	save_figure(figure, path)


def plot_morfr(associations, rankings_path, performance_path, path):
	rankings = pd.read_csv(rankings_path)
	performance = pd.read_csv(performance_path)
	associations['Gene'] = associations['GENES'].fillna('').str.split(';')
	gene_associations = associations.explode('Gene')
	gene_tier = (
		gene_associations.sort_values('q')
		.drop_duplicates('Gene')
		.set_index('Gene')['Tier']
	)
	associated_targets = gene_associations.groupby('Gene')['Phenotype'].agg(lambda values: set(values))
	ranked_genes = rankings[rankings['Gene'].isin(associated_targets.index)].copy()
	ranked_genes['Matching_Target'] = ranked_genes.apply(
		lambda row: row['Target'] in associated_targets[row['Gene']], axis=1
	)
	ranked_genes['Tier'] = ranked_genes['Gene'].map(gene_tier)
	ranked_genes = ranked_genes.sort_values(['MORFR_Rank', 'Gene'])
	ranked_genes['Label'] = ranked_genes.apply(
		lambda row: f"{row['Gene']} / {row['Target']}" + ('' if row['Matching_Target'] else ' (cross-target)'),
		axis=1,
	)
	figure, (rank_axis, score_axis) = plt.subplots(
		1, 2, figsize=(11.5, max(4.8, len(ranked_genes) * 0.36)),
		gridspec_kw={'width_ratios': [1.05, 1.2]},
	)
	for y, (_, row) in enumerate(ranked_genes.iterrows()):
		if row['Matching_Target']:
			rank_axis.scatter(
				row['MORFR_Rank'], y, s=55, color=TIER_COLORS[row['Tier']],
				marker='o', edgecolor='white', linewidth=0.7, zorder=3,
			)
		else:
			rank_axis.scatter(
				row['MORFR_Rank'], y, s=55, facecolor='white',
				edgecolor=TIER_COLORS[row['Tier']], linewidth=1.5,
				marker='s', zorder=3,
			)
		rank_axis.text(row['MORFR_Rank'] + 18, y, str(int(row['MORFR_Rank'])), va='center', fontsize=7.5)
	rank_axis.set_yticks(np.arange(len(ranked_genes)))
	rank_axis.set_yticklabels(ranked_genes['Label'], fontsize=8)
	rank_axis.invert_yaxis()
	rank_axis.set_xlim(0, 1110)
	rank_axis.set_xlabel('MORFR importance rank (1 = highest)')
	rank_axis.set_title('Variant-associated genes in saved top 1,000', loc='left', fontweight='bold')
	rank_axis.axvline(1000, color='#777777', linestyle=':', linewidth=0.8)
	rank_axis.grid(axis='x', color='#eeeeee', linewidth=0.6)
	rank_axis.spines[['top', 'right']].set_visible(False)

	performance = performance.copy()
	performance['Target'] = pd.Categorical(
		performance['Target'],
		categories=['Age', 'SjD Diagnosis', 'Focus Score', 'SSA', 'SSB', 'ANA', 'RF'],
		ordered=True,
	)
	performance = performance.sort_values('Target').reset_index(drop=True)
	y = np.arange(len(performance))
	score_axis.errorbar(
		performance['Outer_holdout_corr2_mean'], y,
		xerr=performance['Outer_holdout_corr2_sd'],
		fmt='o', color='#315f70', ecolor='#7897a0',
		markersize=5.5, capsize=2, linewidth=1.2,
	)
	score_axis.set_yticks(y)
	score_axis.set_yticklabels(performance['Target'].astype(str), fontsize=8)
	score_axis.invert_yaxis()
	score_axis.set_xlim(0, max(0.65, (performance['Outer_holdout_corr2_mean'] + performance['Outer_holdout_corr2_sd']).max() + 0.06))
	score_axis.set_xlabel('Outer-holdout correlation squared (mean +/- SD)')
	score_axis.set_title('Overall target prediction', loc='left', fontweight='bold')
	score_axis.grid(axis='x', color='#eeeeee', linewidth=0.6)
	score_axis.spines[['top', 'right']].set_visible(False)

	variant_tiers_in_figure = set(ranked_genes['Tier'])
	legend_handles = [
		Line2D([0], [0], marker='o', linestyle='none', color='none',
		       markerfacecolor=TIER_COLORS[tier], markeredgecolor='white', markersize=6, label=label)
		for tier, label in [
			('Strong (q<0.05)', 'Variant q<0.05'),
			('Suggestive (0.05<=q<0.10)', 'Variant 0.05<=q<0.10'),
			('Exploratory (0.10<=q<0.25)', 'Variant 0.10<=q<0.25'),
		]
		if tier in variant_tiers_in_figure
	]
	legend_handles.extend([
		Line2D([0], [0], marker='o', linestyle='none', color='#555555',
		       markerfacecolor='#555555', markersize=6, label='Matching phenotype'),
		Line2D([0], [0], marker='s', linestyle='none', color='#555555',
		       markerfacecolor='white', markersize=6, label='Cross-target'),
	])
	figure.legend(handles=legend_handles, frameon=False, ncol=5, loc='upper center',
	              bbox_to_anchor=(0.60, 0.995), fontsize=7)
	figure.subplots_adjust(left=0.24, right=0.98, top=0.82, bottom=0.08, wspace=0.42)
	save_figure(figure, path)


def plot_de_volcanoes(de_dir, path):
	frames = {}
	for label, filename in DE_FILES.items():
		data = pd.read_csv(de_dir / filename, index_col=0)
		data = data[['log2FoldChange', 'padj']].apply(pd.to_numeric, errors='coerce').dropna()
		data['neg_log10_padj'] = -np.log10(data['padj'].clip(lower=1e-300))
		data['significant'] = data['padj'] < 0.05
		frames[label] = data

	all_lfc = pd.concat([frame['log2FoldChange'].abs() for frame in frames.values()])
	x_limit = max(4.0, float(np.ceil(all_lfc.quantile(0.999))))
	max_y = max(float(frame['neg_log10_padj'].max()) for frame in frames.values())
	y_limit = max(8.0, float(np.ceil(max_y + 0.5)))
	figure, axes = plt.subplots(2, 3, figsize=(11.2, 7.4), sharex=True, sharey=True)
	axes = axes.ravel()
	colors = {'not_significant': '#c9ced2', 'down': '#3976a8', 'up': '#bf5148'}
	for axis, (label, data) in zip(axes, frames.items()):
		nonsignificant = data[~data['significant']]
		significant = data[data['significant']]
		axis.scatter(
			nonsignificant['log2FoldChange'], nonsignificant['neg_log10_padj'],
			s=5, color=colors['not_significant'], alpha=0.35, linewidths=0, rasterized=True,
		)
		up = significant[significant['log2FoldChange'] > 0]
		down = significant[significant['log2FoldChange'] <= 0]
		axis.scatter(up['log2FoldChange'], up['neg_log10_padj'], s=6, color=colors['up'], alpha=0.58, linewidths=0, rasterized=True)
		axis.scatter(down['log2FoldChange'], down['neg_log10_padj'], s=6, color=colors['down'], alpha=0.58, linewidths=0, rasterized=True)
		axis.axhline(-np.log10(0.05), color='#333333', linestyle='--', linewidth=0.8)
		axis.axvline(0, color='#777777', linewidth=0.55)
		axis.axvline(-1, color='#888888', linestyle=':', linewidth=0.7)
		axis.axvline(1, color='#888888', linestyle=':', linewidth=0.7)
		axis.set_title(f'{label}  |  {len(significant):,} genes at padj < 0.05', fontsize=9, loc='left')
		axis.set_xlim(-x_limit, x_limit)
		axis.set_ylim(0, y_limit)
		axis.spines[['top', 'right']].set_visible(False)
		axis.grid(axis='y', color='#eeeeee', linewidth=0.5)
	for axis in axes[3:5]:
		axis.set_xlabel('DESeq2 log2 fold-change')
	for axis in (axes[0], axes[3]):
		axis.set_ylabel('-log10 per-contrast DESeq2 adjusted p-value')
	axes[-1].axis('off')
	legend = [
		Line2D([0], [0], marker='o', color='none', markerfacecolor=colors['up'], markersize=5, label='Higher in positive group'),
		Line2D([0], [0], marker='o', color='none', markerfacecolor=colors['down'], markersize=5, label='Lower in positive group'),
		Line2D([0], [0], marker='o', color='none', markerfacecolor=colors['not_significant'], markersize=5, label='padj >= 0.05'),
	]
	figure.legend(handles=legend, frameon=False, ncol=3, loc='lower center', bbox_to_anchor=(0.5, 0.015), fontsize=8)
	figure.suptitle('Gene-level differential expression', x=0.08, ha='left', fontsize=12, fontweight='bold', y=0.99)
	figure.tight_layout(rect=(0.04, 0.05, 0.99, 0.95), h_pad=1.7, w_pad=1.5)
	save_figure(figure, path)


def main():
	parser = argparse.ArgumentParser(description='Generate publication-ready DE, variant-tier, and MORFR figures.')
	parser.add_argument('--scan', type=Path, default=Path('coding_exon_sequence_allele_scan_20261002/all_tests.csv.gz'))
	parser.add_argument('--tier-join', type=Path, default=Path('coding_exon_sequence_allele_scan_20261002/morfr_tier_evidence.csv'))
	parser.add_argument('--morfr-dir', type=Path, default=Path('coding_exon_snp_scan_20261002/morfr_rankings_scored_20261002'))
	parser.add_argument('--de-dir', type=Path, default=Path('de'))
	parser.add_argument('--figure-dir', type=Path, default=Path('figures'))
	parser.add_argument('--overwrite', action='store_true', help='Replace the three named figure outputs in the selected directory.')
	args = parser.parse_args()
	for path in [args.scan, args.tier_join, args.morfr_dir / 'morfr_top1000_by_target.csv', args.morfr_dir / 'morfr_target_performance.csv']:
		if not path.is_file():
			parser.error(f'Input file does not exist: {path}')
	for filename in DE_FILES.values():
		if not (args.de_dir / filename).is_file():
			parser.error(f'DE result does not exist: {args.de_dir / filename}')
	args.figure_dir.mkdir(parents=True, exist_ok=True)
	stems = [
		args.figure_dir / 'variant_association_tiers',
		args.figure_dir / 'morfr_variant_overlap',
		args.figure_dir / 'gene_DESeq2_volcano_panels',
	]
	if not args.overwrite and any(stem.with_suffix(extension).exists() for stem in stems for extension in ('.pdf', '.png')):
		parser.error('At least one publication figure already exists; choose another --figure-dir to avoid overwrite')

	matplotlib.rcParams.update({
		'font.family': 'DejaVu Sans',
		'font.size': 8,
		'axes.labelsize': 8,
		'axes.titlesize': 9,
		'pdf.fonttype': 42,
		'ps.fonttype': 42,
		'axes.linewidth': 0.8,
		'savefig.facecolor': 'white',
	})
	scan = pd.read_csv(args.scan)
	q_columns = [f'{phenotype}_Q_BH_ALL_TESTS' for phenotype in PHENOTYPES]
	associations = load_variant_associations(args.scan)
	plot_variant_tiers(associations, stems[0])
	plot_morfr(
		associations,
		args.morfr_dir / 'morfr_top1000_by_target.csv',
		args.morfr_dir / 'morfr_target_performance.csv',
		stems[1],
	)
	plot_de_volcanoes(args.de_dir, stems[2])
	for stem in stems:
		print(f'Wrote {stem.with_suffix(".pdf")} and {stem.with_suffix(".png")}')


if __name__ == '__main__':
	main()