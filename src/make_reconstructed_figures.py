import argparse
from math import ceil
from pathlib import Path
from textwrap import fill

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Patch, Rectangle
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.decomposition import PCA
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.obo_parser import GODag
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from goatools.test_data.genes_NCBI_9606_All import GENEID2NT

from scan_exon_sequence_alleles import PHENOTYPES, parse_passing_genotype


DE_FILES = {
	'SjD Diagnosis': 'gene.disease.csv',
	'SSA': 'gene.ssa.csv',
	'SSB': 'gene.ssb.csv',
	'ANA': 'gene.ana.csv',
	'RF': 'gene.rf.csv',
}
CLINICAL_TRAITS = [
	'Focus Score', 'Schirmer', 'Tear Film Breakup Time',
	'Rose Bengal/van Bijsterveld', 'Salivary Flow', 'C3', 'C4', 'IgG',
	'Lactate Dehydrogenase', 'Gamma Fraction', 'B2M Fraction',
]
COHORT_COLORS = {'BRA': '#397a78', 'MSG': '#c36f4d'}
OUTCOME_COLORS = {'Negative': '#54728f', 'Positive': '#c36f4d'}
VARIANT_TIER_COLORS = {
	'<0.05': '#b4474b',
	'<0.1': '#d28b19',
	'<0.25': '#3e7884',
}


def save_figure(figure, output_directory, stem):
	figure.savefig(output_directory / f'{stem}.pdf', bbox_inches='tight')
	figure.savefig(output_directory / f'{stem}.png', dpi=300, bbox_inches='tight')
	plt.close(figure)


def q_tier(value):
	if value < 0.05:
		return '<0.05'
	if value < 0.10:
		return '<0.1'
	return '<0.25'


def load_variant_hits(scan_path):
	scan = pd.read_csv(scan_path)
	q_columns = {phenotype: f'{phenotype}_Q_BH_ALL_TESTS' for phenotype in PHENOTYPES}
	missing = [column for column in q_columns.values() if column not in scan]
	if missing:
		raise ValueError('Variant scan is missing pooled q columns: ' + ', '.join(missing))
	minimum_q = scan[list(q_columns.values())].min(axis=1)
	hits = scan[minimum_q < 0.25].copy()
	hits['Min_Q'] = minimum_q[minimum_q < 0.25]
	hits['Lead_Phenotype'] = hits[list(q_columns.values())].idxmin(axis=1).str.replace(
		'_Q_BH_ALL_TESTS', '', regex=False
	)
	hits['Associated_Phenotypes'] = hits.apply(
		lambda row: ';'.join(
			phenotype for phenotype, column in q_columns.items()
			if pd.notna(row[column]) and row[column] < 0.25
		),
		axis=1,
	)
	return hits, q_columns


def extract_candidate_genotypes(hits, metadata, vcf_directory, min_dp=10, min_gq=20):
	sample_ids = metadata.index.astype(str).tolist()
	candidate_by_chrom = {}
	for _, hit in hits.iterrows():
		candidate_by_chrom.setdefault(str(hit['CHROM']), {}).setdefault(int(hit['POS']), []).append(
			(str(hit['REF']), str(hit['ALT']))
		)
	genotypes = {}
	for chrom, positions in candidate_by_chrom.items():
		vcf_path = vcf_directory / f'hs1.{chrom}.vcf'
		if not vcf_path.is_file():
			raise FileNotFoundError(f'Missing VCF: {vcf_path}')
			
		column_index = None
		found = set()
		with vcf_path.open(encoding='utf-8', errors='replace') as vcf:
			for line in vcf:
				if line.startswith('##'):
					continue
				if line.startswith('#CHROM\t'):
					header = line.rstrip('\n').split('\t')
					column_index = {sample: index + 9 for index, sample in enumerate(header[9:])}
					absent = [sample for sample in sample_ids if sample not in column_index]
					if absent:
						raise ValueError(f'{vcf_path.name} lacks metadata samples: {absent[:5]}')
					continue
				if line.startswith('#'):
					continue
				coordinate = line.split('\t', 3)
				position = int(coordinate[1])
				if position not in positions:
					continue
				fields = line.rstrip('\r\n').split('\t')
				ref = fields[3]
				alts = fields[4].split(',')
				format_indices = {key: index for index, key in enumerate(fields[8].split(':'))}
				for candidate_ref, candidate_alt in positions[position]:
					key = (chrom, position, candidate_ref, candidate_alt)
					if key in found or ref != candidate_ref or candidate_alt not in alts:
						continue
					found.add(key)
					allele_index = alts.index(candidate_alt) + 1
					classes = {}
					for sample in sample_ids:
						call = fields[column_index[sample]].split(':')
						parsed = parse_passing_genotype(call, format_indices, min_dp, min_gq)
						group = None
						if parsed is not None:
							alleles, _, _ = parsed
							if allele_index in alleles:
								group = 'ALT carrier'
							elif all(allele == 0 for allele in alleles):
								group = 'Homozygous reference'
						classes[sample] = group
					genotypes[key] = classes
		missing = [
			(chrom, position, ref, alt)
			for position, alleles in positions.items()
			for ref, alt in alleles
			if (chrom, position, ref, alt) not in found
		]
		if missing:
			raise ValueError(f'Could not find {len(missing)} candidate alleles in {vcf_path.name}: {missing[:3]}')
	return genotypes


def genotype_expression_frame(hits, genotype_calls, counts_path, metadata):
	counts = pd.read_csv(counts_path, index_col=0)
	counts.columns = counts.columns.astype(str)
	metadata.index = metadata.index.astype(str)
	if set(counts.columns) != set(metadata.index):
		raise ValueError('Gene-count sample IDs do not exactly match metadata IDs')
	counts = counts.loc[:, metadata.index]
	library_sizes = counts.sum(axis=0).replace(0, np.nan)
	log_cpm = np.log2(counts.divide(library_sizes, axis='columns') * 1_000_000 + 1)
	rows = []
	missing_genes = set()
	for _, hit in hits.iterrows():
		key = (str(hit['CHROM']), int(hit['POS']), str(hit['REF']), str(hit['ALT']))
		groups = genotype_calls[key]
		phenotype = hit['Phenotype']
		for gene in str(hit['GENES']).split(';'):
			if gene not in log_cpm.index:
				missing_genes.add(gene)
				continue
			for sample in metadata.index:
				group = groups.get(sample)
				status = str(metadata.at[sample, phenotype])
				if group is None or status not in OUTCOME_COLORS:
					continue
				rows.append({
					'CHROM': hit['CHROM'], 'POS': hit['POS'], 'REF': hit['REF'], 'ALT': hit['ALT'],
					'VARIANT_TYPE': hit['VARIANT_TYPE'], 'INDEL_LENGTH': hit['INDEL_LENGTH'],
					'Gene': gene, 'Phenotype': phenotype,
					'Variant_Q': hit['Variant_Q'], 'Clinical_Status': status,
					'Sample': sample, 'Cohort': metadata.at[sample, 'Cohort'],
					'Genotype_Group': group, 'Log2_CPM': float(log_cpm.at[gene, sample]),
				})
	return pd.DataFrame(rows), sorted(missing_genes)


def plot_variant_expression(expression, missing_genes, output_directory):
	variant_keys = ['CHROM', 'POS', 'REF', 'ALT', 'Gene', 'Phenotype']
	panels = expression[variant_keys + ['Variant_Q', 'VARIANT_TYPE', 'INDEL_LENGTH']].drop_duplicates()
	panels = panels.sort_values(['Variant_Q', 'Gene', 'Phenotype', 'CHROM', 'POS']).reset_index(drop=True)
	column_count = 3
	row_count = ceil(len(panels) / column_count)
	figure, axes = plt.subplots(
		row_count, column_count,
		figsize=(12.2, max(7.5, row_count * 2.65)),
		squeeze=False,
	)
	rng = np.random.default_rng(20261002)
	base_positions = [0, 1, 3, 4]
	positions = [position + offset for position in base_positions for offset in (-0.18, 0.18)]
	cohort_order = ['BRA', 'MSG']
	group_order = ['Homozygous reference', 'ALT carrier']
	status_order = ['Negative', 'Positive']
	for panel_index, panel in panels.iterrows():
		axis = axes.flat[panel_index]
		key_mask = np.ones(len(expression), dtype=bool)
		for column in variant_keys:
			key_mask &= expression[column].astype(str).to_numpy() == str(panel[column])
		panel_data = expression.loc[key_mask]
		values_by_group = []
		box_colors = []
		tick_labels = []
		for cohort in cohort_order:
			for group in group_order:
				counts_by_status = []
				for status in status_order:
					group_data = panel_data[
						(panel_data['Cohort'] == cohort)
						& (panel_data['Genotype_Group'] == group)
						& (panel_data['Clinical_Status'] == status)
					]
					values_by_group.append(group_data['Log2_CPM'].to_numpy(dtype=float))
					counts_by_status.append(len(group_data))
					box_colors.append(OUTCOME_COLORS[status])
				short_group = 'ref' if group == 'Homozygous reference' else 'ALT'
				tick_labels.append(f'{cohort} {short_group}\nn-/n+ {counts_by_status[0]}/{counts_by_status[1]}')
		box = axis.boxplot(
			values_by_group,
			positions=positions,
			widths=0.22,
			patch_artist=True,
			showfliers=False,
			medianprops={'color': '#202020', 'linewidth': 1.0},
			whiskerprops={'color': '#555555', 'linewidth': 0.7},
			capprops={'color': '#555555', 'linewidth': 0.7},
		)
		for patch, color in zip(box['boxes'], box_colors):
			patch.set_facecolor(color)
			patch.set_alpha(0.35)
		for position, values, color in zip(positions, values_by_group, box_colors):
			if len(values):
				jitter = rng.uniform(-0.07, 0.07, size=len(values))
				axis.scatter(position + jitter, values, s=8, color=color, alpha=0.42, linewidths=0, zorder=3)
		axis.set_xticks(base_positions, tick_labels, fontsize=6.4)
		axis.set_title(
			f"{panel['Gene']} | {panel['Phenotype']} | q={panel['Variant_Q']:.3g}\n"
			f"{panel['CHROM']}:{panel['POS']} {panel['VARIANT_TYPE']} {int(panel['INDEL_LENGTH']):+d} bp",
			fontsize=7.2,
			loc='left',
		)
		axis.set_ylabel('log2(CPM + 1)', fontsize=7)
		axis.grid(axis='y', color='#eeeeee', linewidth=0.5)
		axis.spines[['top', 'right']].set_visible(False)
	for axis in axes.flat[len(panels):]:
		axis.axis('off')
	legend = [
		Patch(facecolor=color, edgecolor='none', alpha=0.55, label=status)
		for status, color in OUTCOME_COLORS.items()
	]
	figure.legend(handles=legend, frameon=False, ncol=2, loc='lower center', bbox_to_anchor=(0.5, 0.01), fontsize=8)
	figure.subplots_adjust(left=0.08, right=0.98, top=0.96, bottom=0.12, hspace=0.65, wspace=0.32)
	save_figure(figure, output_directory, 'variant_expression_by_genotype')


def load_de_tables(de_directory):
	results = {}
	for phenotype, filename in DE_FILES.items():
		data = pd.read_csv(de_directory / filename, index_col=0)
		results[phenotype] = data
	return results


def plot_evidence_matrix(hits, q_columns, de_tables, output_directory):
	gene_associations = hits.assign(Gene=hits['GENES'].fillna('').str.split(';')).explode('Gene')
	gene_associations = gene_associations[gene_associations['Gene'] != ''].copy()
	gene_associations['Row_Label'] = gene_associations.apply(
		lambda row: f"{row['Gene']} | {row['CHROM']}:{row['POS']} | {row['VARIANT_TYPE']} {int(row['INDEL_LENGTH']):+d} bp",
		axis=1,
	)
	gene_associations = gene_associations.sort_values(['Min_Q', 'Gene', 'CHROM', 'POS']).reset_index(drop=True)
	row_count = len(gene_associations)
	phenotype_labels = ['SjD', 'SSA', 'SSB', 'ANA', 'RF']
	figure, axes = plt.subplots(
		1, 2,
		figsize=(12.4, max(7.2, row_count * 0.32)),
		sharey=True,
		gridspec_kw={'width_ratios': [1.0, 1.0]},
	)
	q_norm = Normalize(0, 4)
	q_map = plt.get_cmap('YlOrRd')
	for axis in axes:
		axis.set_xlim(0, len(PHENOTYPES))
		axis.set_xticks(np.arange(len(PHENOTYPES)) + 0.5, phenotype_labels, fontsize=8)
		axis.set_yticks(np.arange(row_count) + 0.5, gene_associations['Row_Label'], fontsize=7)
		axis.set_ylim(row_count, 0)
		axis.set_xticks(np.arange(len(PHENOTYPES) + 1), minor=True)
		axis.set_yticks(np.arange(row_count + 1), minor=True)
		axis.grid(which='minor', color='white', linewidth=1.2)
		axis.tick_params(which='minor', bottom=False, left=False)
		axis.spines[['top', 'right', 'left', 'bottom']].set_visible(False)

	for row_index, (_, association) in enumerate(gene_associations.iterrows()):
		gene = association['Gene']
		for phenotype_index, phenotype in enumerate(PHENOTYPES):
			qvalue = association[q_columns[phenotype]]
			if pd.isna(qvalue):
				axes[0].add_patch(Rectangle((phenotype_index, row_index), 1, 1, color='#ececec'))
				axes[0].text(phenotype_index + 0.5, row_index + 0.5, 'NA', ha='center', va='center', fontsize=6, color='#555555')
				continue
			transformed_q = min(-np.log10(max(float(qvalue), 1e-300)), 4)
			axes[0].add_patch(Rectangle((phenotype_index, row_index), 1, 1, color=q_map(q_norm(transformed_q))))
			text_color = 'white' if transformed_q > 2.2 else '#333333'
			axes[0].text(phenotype_index + 0.5, row_index + 0.5, f'{qvalue:.2g}', ha='center', va='center', fontsize=6.2, color=text_color)

	for row_index, association in gene_associations.iterrows():
		gene = association['Gene']
		for phenotype_index, phenotype in enumerate(PHENOTYPES):
			de = de_tables[phenotype]
			if gene not in de.index:
				axes[1].add_patch(Rectangle((phenotype_index, row_index), 1, 1, color='#ececec'))
				axes[1].text(phenotype_index + 0.5, row_index + 0.5, 'NA', ha='center', va='center', fontsize=6, color='#555555')
				continue
			padj = pd.to_numeric(pd.Series([de.at[gene, 'padj']]), errors='coerce').iloc[0]
			lfc = pd.to_numeric(pd.Series([de.at[gene, 'log2FoldChange']]), errors='coerce').iloc[0]
			if pd.isna(lfc):
				axes[1].add_patch(Rectangle((phenotype_index, row_index), 1, 1, color='#ececec'))
				axes[1].text(phenotype_index + 0.5, row_index + 0.5, 'NA', ha='center', va='center', fontsize=6, color='#555555')
				continue
			color = plt.get_cmap('RdBu_r')(Normalize(-2, 2)(np.clip(lfc, -2, 2)))
			axes[1].add_patch(Rectangle((phenotype_index, row_index), 1, 1, color=color))
			label = f'{lfc:+.1f}' + ('*' if pd.notna(padj) and padj < 0.05 else '')
			axes[1].text(phenotype_index + 0.5, row_index + 0.5, label, ha='center', va='center', fontsize=6.2, color='white' if abs(lfc) > 1 else '#333333')

	axes[0].set_title('Pooled CMH q (all phenotypes)', loc='left', fontsize=8.5, fontweight='bold', pad=7)
	axes[1].set_title('DESeq2 log2FC (* padj < 0.05)', loc='left', fontsize=8.5, fontweight='bold', pad=7)
	axes[0].set_ylabel('Variant / overlapping gene')
	figure.colorbar(
		plt.cm.ScalarMappable(norm=q_norm, cmap=q_map),
		ax=axes[0], fraction=0.035, pad=0.025, label='-log10 pooled q (capped at 4)',
	)
	figure.colorbar(
		plt.cm.ScalarMappable(norm=Normalize(-2, 2), cmap='RdBu_r'),
		ax=axes[1], fraction=0.035, pad=0.025, label='log2 fold-change (clipped at +/-2)',
	)
	figure.suptitle('Variant and expression evidence by phenotype', x=0.04, ha='left', y=0.99, fontsize=10, fontweight='bold')
	figure.subplots_adjust(left=0.28, right=0.97, top=0.94, bottom=0.06, wspace=0.46)
	save_figure(figure, output_directory, 'matched_gene_evidence_matrix')


def plot_de_counts(de_tables, output_directory):
	labels = []
	counts = []
	totals = []
	for phenotype in PHENOTYPES:
		data = de_tables[phenotype]
		adjusted = pd.to_numeric(data['padj'], errors='coerce')
		labels.append(phenotype)
		counts.append(int((adjusted < 0.05).sum()))
		totals.append(int(adjusted.notna().sum()))
	figure, axis = plt.subplots(figsize=(7.4, 4.4))
	positions = np.arange(len(PHENOTYPES))
	axis.barh(positions, counts, color='#397a78', height=0.67)
	for position, count, total in zip(positions, counts, totals):
		axis.text(count + max(counts) * 0.012, position, f'{count:,} / {total:,}', va='center', fontsize=8)
	axis.set_yticks(positions, labels)
	axis.invert_yaxis()
	axis.set_xlabel('Genes with per-contrast DESeq2 padj < 0.05')
	axis.set_title('Differential-expression yield by phenotype', loc='left', fontweight='bold')
	axis.set_xlim(0, max(counts) * 1.35)
	axis.grid(axis='x', color='#eeeeee', linewidth=0.6)
	axis.spines[['top', 'right']].set_visible(False)
	figure.tight_layout()
	save_figure(figure, output_directory, 'de_gene_counts_per_contrast')


def expression_log_cpm(counts):
	library_sizes = counts.sum(axis=0).replace(0, np.nan)
	return np.log2(counts.divide(library_sizes, axis='columns') * 1_000_000 + 1)


def plot_expression_pca(counts_path, metadata, output_directory):
	counts = pd.read_csv(counts_path, index_col=0)
	counts.columns = counts.columns.astype(str)
	metadata = metadata.copy()
	metadata.index = metadata.index.astype(str)
	log_cpm = expression_log_cpm(counts.loc[:, metadata.index])
	variable_genes = log_cpm.loc[:, log_cpm.nunique() > 1]
	pca = PCA(n_components=2)
	scores = pca.fit_transform(variable_genes.T)
	plot_data = pd.DataFrame(scores, index=variable_genes.columns, columns=['PC1', 'PC2']).join(metadata)
	variance = pca.explained_variance_ratio_ * 100
	figure, axes = plt.subplots(1, 2, figsize=(9.2, 4.2), sharex=True, sharey=True)
	views = [
		('Cohort', {'BRA': COHORT_COLORS['BRA'], 'MSG': COHORT_COLORS['MSG']}),
		('SjD Diagnosis', {'Negative': '#3976a8', 'Positive': '#c36f4d'}),
	]
	for axis, (column, palette) in zip(axes, views):
		for label, color in palette.items():
			group = plot_data[plot_data[column] == label]
			axis.scatter(group['PC1'], group['PC2'], s=22, color=color, alpha=0.78, linewidths=0.35, edgecolors='white', label=label)
		axis.set_xlabel(f'PC1 ({variance[0]:.1f}%)')
		axis.set_title(column, loc='left', fontweight='bold')
		axis.grid(color='#eeeeee', linewidth=0.5)
		axis.spines[['top', 'right']].set_visible(False)
		axis.legend(frameon=False, fontsize=8)
	axes[0].set_ylabel(f'PC2 ({variance[1]:.1f}%)')
	figure.suptitle('Expression PCA after ComBat-seq adjustment', x=0.06, ha='left', y=1.02, fontsize=11, fontweight='bold')
	figure.tight_layout()
	save_figure(figure, output_directory, 'rna_expression_pca_cohort_diagnosis')


def plot_clinical_measurement_qc(metadata, output_directory):
	values = metadata[CLINICAL_TRAITS].apply(pd.to_numeric, errors='coerce')
	cohorts = metadata['Cohort'].astype(str)
	cohort_sizes = cohorts.value_counts()
	available = pd.DataFrame(index=CLINICAL_TRAITS, columns=['BRA', 'MSG'], dtype=float)
	for cohort in ['BRA', 'MSG']:
		cohort_mask = cohorts == cohort
		available[cohort] = values.loc[cohort_mask].notna().sum(axis=0) / cohort_sizes[cohort]

	rank_values = pd.DataFrame(index=values.index, columns=CLINICAL_TRAITS, dtype=float)
	for cohort in ['BRA', 'MSG']:
		cohort_mask = cohorts == cohort
		rank_values.loc[cohort_mask] = values.loc[cohort_mask].rank(pct=True)
	correlation = rank_values.corr(min_periods=15)

	figure, (coverage_axis, correlation_axis) = plt.subplots(
		1, 2, figsize=(14.0, 7.2), gridspec_kw={'width_ratios': [0.85, 1.4]},
	)
	positions = np.arange(len(CLINICAL_TRAITS))
	bar_height = 0.34
	for offset, cohort in [(-bar_height / 2, 'BRA'), (bar_height / 2, 'MSG')]:
		coverage_axis.barh(
			positions + offset,
			available[cohort].to_numpy(),
			height=bar_height,
			color=COHORT_COLORS[cohort],
			label=cohort,
		)
	coverage_axis.set_yticks(positions, CLINICAL_TRAITS, fontsize=7.5)
	coverage_axis.invert_yaxis()
	coverage_axis.set_xlim(0, 1.08)
	coverage_axis.set_xlabel('Fraction with measurement')
	coverage_axis.set_title('Measurement coverage', loc='left', fontweight='bold')
	coverage_axis.legend(frameon=False, fontsize=8)
	coverage_axis.grid(axis='x', color='#eeeeee', linewidth=0.5)
	coverage_axis.spines[['top', 'right']].set_visible(False)

	image = correlation_axis.imshow(correlation.to_numpy(dtype=float), cmap='RdBu_r', vmin=-1, vmax=1)
	correlation_axis.set_xticks(positions, [trait.replace('Tear Film Breakup Time', 'TBUT').replace('Rose Bengal/van Bijsterveld', 'Rose Bengal').replace('Lactate Dehydrogenase', 'LDH') for trait in CLINICAL_TRAITS], rotation=55, ha='right', fontsize=7)
	correlation_axis.set_yticks(positions, [trait.replace('Tear Film Breakup Time', 'TBUT').replace('Rose Bengal/van Bijsterveld', 'Rose Bengal').replace('Lactate Dehydrogenase', 'LDH') for trait in CLINICAL_TRAITS], fontsize=7)
	correlation_axis.set_title('Within-cohort rank correlation', loc='left', fontweight='bold')
	for row in range(len(CLINICAL_TRAITS)):
		for column in range(len(CLINICAL_TRAITS)):
			value = correlation.iat[row, column]
			if pd.notna(value):
				correlation_axis.text(column, row, f'{value:.1f}', ha='center', va='center', fontsize=6, color='white' if abs(value) > 0.55 else '#333333')
	figure.colorbar(image, ax=correlation_axis, fraction=0.046, pad=0.04, label='Spearman rho')
	figure.suptitle('Clinical data coverage and covariation', x=0.04, ha='left', y=0.99, fontsize=11, fontweight='bold')
	figure.subplots_adjust(left=0.20, right=0.98, top=0.91, bottom=0.18, wspace=0.40)
	save_figure(figure, output_directory, 'clinical_measurement_qc')


def run_corrected_go_enrichment(scan_path, counts_path, de_tables, rankings_path, obo_path, gene2go_path):
	gene2go = Gene2GoReader(str(gene2go_path), taxids=[9606]).get_ns2assc()
	go_dag = GODag(str(obo_path))
	counts = pd.read_csv(counts_path, index_col=0)
	cpm = counts.divide(counts.sum(axis=0).replace(0, np.nan), axis='columns') * 1_000_000
	expressed_symbols = set(counts.index[(cpm >= 1).sum(axis=1) >= ceil(counts.shape[1] * 0.10)])
	symbol_to_entrez = {}
	for entrez_id, gene_record in GENEID2NT.items():
		symbol_to_entrez.setdefault(str(gene_record[5]), set()).add(int(entrez_id))
	background = set().union(*(symbol_to_entrez.get(symbol, set()) for symbol in expressed_symbols))
	background.intersection_update(GENEID2NT.keys())
	goea = GOEnrichmentStudyNS(background, gene2go, go_dag, alpha=0.05, methods=['fdr_bh'])

	scan = pd.read_csv(scan_path)
	rankings = pd.read_csv(rankings_path)
	output_rows = []
	for phenotype in PHENOTYPES:
		q_column = f'{phenotype}_Q_BH_ALL_TESTS'
		variant_genes = set()
		for gene_cell in scan.loc[scan[q_column] < 0.25, 'GENES'].dropna():
			variant_genes.update(gene_cell.split(';'))
		de_table = de_tables[phenotype]
		de_adjusted = pd.to_numeric(de_table['padj'], errors='coerce')
		de_genes = set(de_table.index[de_adjusted < 0.05])
		morfr_genes = set(
			rankings.loc[rankings['Target'] == phenotype]
			.sort_values('MORFR_Rank').head(1000)['Gene'].astype(str)
		)
		gene_sets = {
			'DESeq2': de_genes,
			'Variant q<0.25': variant_genes,
			'MORFR top 1,000': morfr_genes,
		}
		for evidence_source, genes in gene_sets.items():
			study_ids = set().union(*(symbol_to_entrez.get(gene, set()) for gene in genes))
			study_ids.intersection_update(background)
			if len(study_ids) < 2:
				continue
			results = goea.run_study(study_ids, prt=None)
			for result in results:
				if result.p_fdr_bh >= 0.05:
					continue
				output_rows.append({
					'Phenotype': phenotype,
					'Evidence_Source': evidence_source,
					'GO_ID': result.GO,
					'Namespace': result.NS,
					'Term': result.name,
					'Study_Count': result.study_count,
					'Study_Size': result.study_n,
					'Background_Count': result.pop_count,
					'Background_Size': result.pop_n,
					'Q_BH_Within_Set': result.p_fdr_bh,
				})
	return pd.DataFrame(output_rows)


def plot_go_evidence(enrichment, output_directory):
	source_order = ['DESeq2', 'MORFR top 1,000']
	figure, axes = plt.subplots(len(PHENOTYPES), 1, figsize=(12.5, 13.0), squeeze=False)
	axes = axes.ravel()
	color_map = plt.get_cmap('viridis')
	color_norm = Normalize(0.5, 2.5)
	for axis, phenotype in zip(axes, PHENOTYPES):
		panel = enrichment[enrichment['Phenotype'] == phenotype]
		morfr_terms = panel[panel['Evidence_Source'] == 'MORFR top 1,000']
		term_order = morfr_terms.sort_values('Q_BH_Within_Set').drop_duplicates('GO_ID').head(6)
		de_terms = panel[panel['Evidence_Source'] == 'DESeq2']
		shared_terms = set(morfr_terms['GO_ID']) & set(de_terms['GO_ID'])
		for row_index, term in enumerate(term_order.itertuples()):
			for source_index, source in enumerate(source_order):
				match = panel[(panel['GO_ID'] == term.GO_ID) & (panel['Evidence_Source'] == source)]
				if match.empty:
					continue
				result = match.iloc[0]
				gene_ratio = result['Study_Count'] / result['Study_Size'] if result['Study_Size'] else 0
				qvalue = max(float(result['Q_BH_Within_Set']), 1e-300)
				color_value = np.log10(1 + -np.log10(qvalue))
				axis.scatter(
					source_index, row_index,
					s=32 + min(110, 700 * gene_ratio),
					c=[color_value], cmap=color_map, norm=color_norm,
					edgecolor='white', linewidth=0.6,
				)
		axis.set_yticks(
			np.arange(len(term_order)),
			[term.Term + f' [{term.Namespace}]' for term in term_order.itertuples()],
			fontsize=7,
		)
		axis.invert_yaxis()
		axis.set_xticks([0, 1], ['DESeq2', 'MORFR'], fontsize=7.5)
		axis.set_xlim(-0.5, 1.5)
		axis.set_title(
			f'{phenotype} | shared terms {len(shared_terms)}/{morfr_terms.GO_ID.nunique()} MORFR terms',
			loc='left', fontsize=8.5, fontweight='bold', pad=3,
		)
		axis.grid(axis='y', color='#eeeeee', linewidth=0.4)
		axis.spines[['top', 'right', 'left', 'bottom']].set_visible(False)
		axis.tick_params(axis='y', length=0)
	colorbar = figure.colorbar(
		plt.cm.ScalarMappable(norm=color_norm, cmap=color_map),
		ax=axes, fraction=0.025, pad=0.02,
	)
	colorbar.set_label('log10(1 + -log10(q))')
	figure.suptitle('Shared GO terms in phenotype DE and MORFR lists', x=0.04, ha='left', y=0.99, fontsize=11, fontweight='bold')
	figure.text(0.04, 0.955, 'Rows are ordered by MORFR q and shared across DE/MORFR columns. No variant-overlap set had GO q<0.05.', fontsize=8, color='#333333')
	figure.subplots_adjust(left=0.35, right=0.92, top=0.92, bottom=0.05, hspace=0.44)
	save_figure(figure, output_directory, 'go_enrichment_by_evidence')


def main():
	parser = argparse.ArgumentParser(description='Reconstruct legacy res9, Venn, and cohort-QC figure stories from corrected analyses.')
	parser.add_argument('--scan', type=Path, default=Path('coding_exon_sequence_allele_scan_20261002/all_tests.csv.gz'))
	parser.add_argument('--counts', type=Path, default=Path('all_gene.csv'))
	parser.add_argument('--metadata', type=Path, default=Path('meta.csv'))
	parser.add_argument('--de-dir', type=Path, default=Path('de'))
	parser.add_argument('--morfr-dir', type=Path, default=Path('coding_exon_snp_scan_20261002/morfr_rankings_scored_20261002'))
	parser.add_argument('--vcf-dir', type=Path, default=Path(r'C:\Users\thoma\Dropbox\HelixBiowulf\2021-02-03_GATK'))
	parser.add_argument('--go-obo', type=Path, default=Path(r'..\go-basic.obo'))
	parser.add_argument('--gene2go', type=Path, default=Path(r'..\gene2go'))
	parser.add_argument('--go-table', type=Path, default=Path('coding_exon_sequence_allele_scan_20261002/go_enrichment_by_evidence.csv'))
	parser.add_argument('--figure-dir', type=Path, default=Path('figures/publication_20261002'))
	parser.add_argument('--overwrite', action='store_true')
	args = parser.parse_args()
	metadata = pd.read_csv(args.metadata, index_col=0)
	for path in [
		args.scan, args.counts, args.metadata,
		args.morfr_dir / 'morfr_top1000_by_target.csv',
		args.go_obo, args.gene2go,
	]:
		if not path.is_file():
			parser.error(f'Input file does not exist: {path}')
	for filename in DE_FILES.values():
		if not (args.de_dir / filename).is_file():
			parser.error(f'DE result does not exist: {args.de_dir / filename}')
	output_stems = [
		'variant_expression_by_genotype', 'matched_gene_evidence_matrix',
		'de_gene_counts_per_contrast', 'rna_expression_pca_cohort_diagnosis',
		'clinical_measurement_qc', 'go_enrichment_by_evidence',
	]
	if not args.overwrite and any(
		(args.figure_dir / f'{stem}.{extension}').exists()
		for stem in output_stems for extension in ['pdf', 'png']
	):
		parser.error('At least one reconstructed figure exists; pass --overwrite to replace the set')
	args.figure_dir.mkdir(parents=True, exist_ok=True)
	matplotlib.rcParams.update({
		'font.family': 'DejaVu Sans', 'font.size': 8,
		'axes.labelsize': 8, 'axes.titlesize': 9,
		'pdf.fonttype': 42, 'ps.fonttype': 42,
		'axes.linewidth': 0.8, 'savefig.facecolor': 'white',
	})
	hits, q_columns = load_variant_hits(args.scan)
	phenotype_hits = hits.melt(
		id_vars=['CHROM', 'POS', 'REF', 'ALT', 'VARIANT_TYPE', 'INDEL_LENGTH', 'GENES'],
		value_vars=list(q_columns.values()),
		var_name='Phenotype_Q',
		value_name='Variant_Q',
	)
	phenotype_hits = phenotype_hits[phenotype_hits['Variant_Q'] < 0.05].copy()
	phenotype_hits['Phenotype'] = phenotype_hits['Phenotype_Q'].str.replace(
		'_Q_BH_ALL_TESTS', '', regex=False
	)
	strong_variants = phenotype_hits.drop_duplicates(['CHROM', 'POS', 'REF', 'ALT'])
	genotypes = extract_candidate_genotypes(strong_variants, metadata, args.vcf_dir)
	expression, missing_genes = genotype_expression_frame(phenotype_hits, genotypes, args.counts, metadata)
	if expression.empty:
		parser.error('No expression observations were available for q-tier candidate genes')
	plot_variant_expression(expression, missing_genes, args.figure_dir)
	de_tables = load_de_tables(args.de_dir)
	plot_evidence_matrix(hits, q_columns, de_tables, args.figure_dir)
	plot_de_counts(de_tables, args.figure_dir)
	plot_expression_pca(args.counts, metadata, args.figure_dir)
	plot_clinical_measurement_qc(metadata, args.figure_dir)
	go_enrichment = run_corrected_go_enrichment(
		args.scan, args.counts, de_tables,
		args.morfr_dir / 'morfr_top1000_by_target.csv',
		args.go_obo, args.gene2go,
	)
	args.go_table.parent.mkdir(parents=True, exist_ok=True)
	go_enrichment.to_csv(args.go_table, index=False)
	plot_go_evidence(go_enrichment, args.figure_dir)
	print(f'Expression comparison: {expression.Gene.nunique()} genes, {expression.Sample.nunique()} samples, {len(missing_genes)} overlapping annotation IDs without gene-count rows')
	print(f'Wrote corrected GO terms: {args.go_table} ({len(go_enrichment)} terms at within-set q<0.05)')
	for stem in output_stems:
		print(f'Wrote {args.figure_dir / (stem + ".pdf")} and {args.figure_dir / (stem + ".png")}')


if __name__ == '__main__':
	main()