import argparse
from collections import defaultdict
from pathlib import Path

import pandas as pd


PHENOTYPES = [
	('SjD Diagnosis', 'disease', 'SjD_Diagnosis'),
	('SSA', 'ssa', 'SSA'),
	('SSB', 'ssb', 'SSB'),
	('ANA', 'ana', 'ANA'),
	('RF', 'rf', 'RF'),
]
SITE_KEY = ['#CHROM', 'POS', 'REF', 'ALT']


def gene_symbol(transcript_name):
	return transcript_name.rsplit('-', 1)[0]


def annotate_exon_overlaps(candidates, bed_path):
	by_chrom = defaultdict(list)
	for index, row in candidates.iterrows():
		by_chrom[str(row['#CHROM'])].append((index, int(row['POS']) - 1))
	genes = {index: set() for index in candidates.index}
	with bed_path.open(encoding='utf-8') as bed:
		for line in bed:
			fields = line.rstrip('\n').split('\t')
			if len(fields) < 9 or not fields[1].isdigit() or fields[7] != 'protein_coding':
				continue
			targets = by_chrom.get(fields[0], [])
			if not targets:
				continue
			transcript_start = int(fields[1])
			block_sizes = [int(value) for value in fields[5].rstrip(',').split(',') if value]
			block_offsets = [int(value) for value in fields[6].rstrip(',').split(',') if value]
			if len(block_sizes) != len(block_offsets):
				raise ValueError(f'Mismatched exon-block arrays for {fields[3]}')
			exons = [
				(transcript_start + offset, transcript_start + offset + size)
				for size, offset in zip(block_sizes, block_offsets)
			]
			name = gene_symbol(fields[3])
			for index, position in targets:
				if any(start <= position < end for start, end in exons):
					genes[index].add(name)
	return genes


def main():
	parser = argparse.ArgumentParser(
		description='Build an exploratory variant, DE, and descriptive expression evidence table.'
	)
	parser.add_argument('--associations', type=Path, default=Path('coding_exon_snp_scan_20261002/all_associations_verified.csv.gz'))
	parser.add_argument('--bed', type=Path, default=Path(r'C:\Users\thoma\Dropbox\HelixBiowulf\2021-02-03_GATK\all_ens.bed'))
	parser.add_argument('--de-dir', type=Path, default=Path('de'))
	parser.add_argument('--var1', type=Path, default=Path('var1.csv'))
	parser.add_argument('--meta', type=Path, default=Path('meta.csv'))
	parser.add_argument('--q-cutoff', type=float, default=0.25)
	parser.add_argument('--out', type=Path, default=Path('coding_exon_snp_scan_20261002/exploratory_evidence_q25.csv'))
	args = parser.parse_args()
	if args.out.exists():
		parser.error(f'Output already exists; choose a new path: {args.out}')
	if not 0 < args.q_cutoff <= 1:
		parser.error('--q-cutoff must be in (0, 1]')
	for path in (args.associations, args.bed, args.var1, args.meta):
		if not path.is_file():
			parser.error(f'Input file does not exist: {path}')

	associations = pd.read_csv(args.associations)
	all_five_q_columns = [f'{phenotype}_Q_CMH_ALL_PHENOTYPES' for phenotype, _, _ in PHENOTYPES]
	missing = [column for column in SITE_KEY + all_five_q_columns if column not in associations]
	if missing:
		parser.error('Association file lacks columns: ' + ', '.join(missing))
	all_five_q = associations[all_five_q_columns].apply(pd.to_numeric, errors='coerce')
	minimum_q = all_five_q.min(axis=1)
	candidates = associations.loc[minimum_q <= args.q_cutoff].copy()
	candidates['Min_Q_All_Phenotypes'] = minimum_q.loc[candidates.index]
	candidates['Evidence_Tier'] = pd.cut(
		candidates['Min_Q_All_Phenotypes'],
		bins=[-float('inf'), 0.05, 0.10, args.q_cutoff, float('inf')],
		labels=['q<=0.05_within_scan', 'q<=0.10_suggestive', 'q<=cutoff_exploratory', ''],
		right=True,
	)
	candidates['Lead_Phenotype'] = all_five_q.loc[candidates.index].idxmin(axis=1).str.replace(
		'_Q_CMH_ALL_PHENOTYPES', '', regex=False
	)

	genes_by_candidate = annotate_exon_overlaps(candidates, args.bed)
	var1 = pd.read_csv(args.var1, index_col=0)
	de_tables = {
		phenotype: pd.read_csv(args.de_dir / f'gene.{de_key}.csv', index_col=0)
		for phenotype, de_key, _ in PHENOTYPES
	}
	metadata = pd.read_csv(args.meta, index_col=0)

	metadata_counts = {}
	for phenotype, _, _ in PHENOTYPES:
		for cohort in sorted(metadata['Cohort'].dropna().unique()):
			values = metadata.loc[metadata['Cohort'] == cohort, phenotype]
			metadata_counts[f'Meta_{phenotype}_{cohort}_N_positive'] = int((values == 'Positive').sum())
			metadata_counts[f'Meta_{phenotype}_{cohort}_N_negative'] = int((values == 'Negative').sum())
			metadata_counts[f'Meta_{phenotype}_{cohort}_N_missing'] = int(values.isna().sum())

	rows = []
	for index, candidate in candidates.iterrows():
		genes = sorted(genes_by_candidate[index])
		if not genes:
			genes = ['']
		for gene in genes:
			row = {
				column: candidate[column]
				for column in [
					*SITE_KEY, 'QUAL', 'FILTER', 'N_PASS', 'PASS_RATE',
					'N_ALT_CARRIERS', 'N_REF', 'ALT_ALLELE_FREQUENCY',
					'MEDIAN_DP', 'MEDIAN_GQ', 'Min_Q_All_Phenotypes',
					'Evidence_Tier', 'Lead_Phenotype',
				]
				if column in candidate
			}
			row['Gene'] = gene
			row['Gene_Overlap'] = 'protein_coding_exon' if gene else 'unannotated'
			for phenotype, _, var1_column in PHENOTYPES:
				for suffix in ['P_FISHER', 'P_CMH', 'Q_FISHER', 'Q_CMH', 'P_MISSING_CMH', 'Q_MISSING_CMH', 'Q_CMH_ALL_PHENOTYPES']:
					column = f'{phenotype}_{suffix}'
					if column in candidate:
						row[column] = candidate[column]
				de = de_tables[phenotype]
				de_match = de.loc[gene] if gene and gene in de.index else None
				for statistic in ['baseMean', 'log2FoldChange', 'pvalue', 'padj']:
					row[f'DE_{phenotype}_{statistic}'] = de_match.get(statistic, pd.NA) if de_match is not None else pd.NA
				row[f'var1_rho2_all_samples_{phenotype}'] = (
					var1.at[gene, var1_column] if gene and gene in var1.index and var1_column in var1.columns else pd.NA
				)
			row.update(metadata_counts)
			row['RF_Importance_Available'] = False
			rows.append(row)

	result = pd.DataFrame(rows)
	args.out.parent.mkdir(parents=True, exist_ok=True)
	result.to_csv(args.out, index=False)
	print(f'Unique candidate sites at q<={args.q_cutoff}: {len(candidates)}')
	print(f'Variant-gene rows written: {len(result)}')
	print('var1 values are descriptive all-sample rho^2, not held-out tree importances; no RF Summary artifacts were found.')
	print('BED annotation is protein-coding transcript exon overlap, not CDS/protein-altering consequence.')
	print(f'Wrote {args.out}')


if __name__ == '__main__':
	main()