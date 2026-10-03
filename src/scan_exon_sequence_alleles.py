import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests


PHENOTYPES = ['SjD Diagnosis', 'SSA', 'SSB', 'ANA', 'RF']
CODING_GENE_TYPES = {'protein_coding'}
CHROMOSOMES = [
	'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
	'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
	'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrM', 'chrX', 'chrY',
]


def number_or_none(value, cast):
	if value in ('', '.', 'NA'):
		return None
	try:
		return cast(value)
	except ValueError:
		return None


def is_sequence_allele(allele):
	return bool(allele) and all(base in 'ACGT' for base in allele)


def variant_type(ref, alt):
	if len(ref) == len(alt) == 1:
		return 'SNV'
	if len(ref) == len(alt):
		return 'MNV'
	return 'insertion' if len(alt) > len(ref) else 'deletion'


def load_coding_exon_intervals(path):
	intervals = {}
	with path.open(encoding='utf-8') as bed:
		for line in bed:
			fields = line.rstrip('\n').split('\t')
			if len(fields) < 9 or not fields[1].isdigit():
				continue
			if fields[7] not in CODING_GENE_TYPES and not fields[7].startswith(('IG_', 'TR_')):
				continue
			transcript = fields[3]
			gene = re.sub(r'-\d+$', '', transcript)
			start = int(fields[1])
			sizes = [int(value) for value in fields[5].rstrip(',').split(',') if value]
			offsets = [int(value) for value in fields[6].rstrip(',').split(',') if value]
			if len(sizes) != len(offsets):
				raise ValueError(f'Mismatched exon-block arrays for {transcript}')
			for size, offset in zip(sizes, offsets):
				intervals.setdefault(fields[0], []).append((start + offset, start + offset + size, gene))

	merged = {}
	for chrom, blocks in intervals.items():
		blocks.sort()
		chrom_blocks = []
		for start, end, gene in blocks:
			if chrom_blocks and start < chrom_blocks[-1][1]:
				previous_start, previous_end, genes = chrom_blocks[-1]
				genes.add(gene)
				chrom_blocks[-1] = (previous_start, max(previous_end, end), genes)
			else:
				chrom_blocks.append((start, end, {gene}))
		starts = np.fromiter((block[0] for block in chrom_blocks), dtype=np.int64)
		merged[chrom] = (starts, chrom_blocks)
	return merged


def overlapping_genes(intervals, chrom, start, end):
	entry = intervals.get(chrom)
	if entry is None:
		return []
	starts, blocks = entry
	stop = int(np.searchsorted(starts, end, side='left'))
	first = max(0, int(np.searchsorted(starts, start, side='left')) - 1)
	genes = set()
	for block_start, block_end, block_genes in blocks[first:stop]:
		if block_start < end and block_end > start:
			genes.update(block_genes)
	return sorted(genes)


def parse_passing_genotype(sample_fields, format_indices, min_dp, min_gq):
	gt_index = format_indices.get('GT')
	if gt_index is None or gt_index >= len(sample_fields):
		return None
	alleles = sample_fields[gt_index].replace('|', '/').split('/')
	if not alleles or not all(allele.isdigit() for allele in alleles):
		return None
	dp_index = format_indices.get('DP')
	gq_index = format_indices.get('GQ')
	depth = number_or_none(sample_fields[dp_index], int) if dp_index is not None and dp_index < len(sample_fields) else None
	quality = number_or_none(sample_fields[gq_index], float) if gq_index is not None and gq_index < len(sample_fields) else None
	if min_dp and (depth is None or depth < min_dp):
		return None
	if min_gq and (quality is None or quality < min_gq):
		return None
	return tuple(int(allele) for allele in alleles), depth, quality


def cmh_pvalue(tables):
	deviation = 0.0
	variance = 0.0
	for table in tables:
		counts = np.asarray(table, dtype=float)
		row_totals = counts.sum(axis=1)
		column_totals = counts.sum(axis=0)
		total = counts.sum()
		if total <= 1 or np.any(row_totals == 0) or np.any(column_totals == 0):
			continue
		deviation += counts[1, 1] - row_totals[1] * column_totals[1] / total
		variance += (
			row_totals[0] * row_totals[1] * column_totals[0] * column_totals[1]
			/ (total * total * (total - 1))
		)
	if variance <= 0:
		return np.nan
	return float(chi2.sf(deviation * deviation / variance, 1))


def allele_tables(carriers, reference, tested, phenotype, cohorts):
	tables = []
	for cohort in np.unique(cohorts):
		in_cohort = tested & (cohorts == cohort) & np.isin(phenotype, ['Negative', 'Positive'])
		negative = in_cohort & (phenotype == 'Negative')
		positive = in_cohort & (phenotype == 'Positive')
		tables.append([
			[int(np.sum(reference & negative)), int(np.sum(reference & positive))],
			[int(np.sum(carriers & negative)), int(np.sum(carriers & positive))],
		])
	return tables


def count_table(tables):
	return [
		sum(table[row][column] for table in tables)
		for row, column in [(0, 0), (0, 1), (1, 0), (1, 1)]
	]


def scan_vcf(vcf_path, sample_ids, cohorts, phenotypes, intervals, args):
	column_index = None
	rows = []
	with vcf_path.open(encoding='utf-8', errors='replace') as vcf:
		for line in vcf:
			if line.startswith('#'):
				if line.startswith('#CHROM\t'):
					header = line.rstrip('\n').split('\t')
					column_index = {sample: index + 9 for index, sample in enumerate(header[9:])}
					absent = [sample for sample in sample_ids if sample not in column_index]
					if absent:
						raise ValueError(f'{vcf_path.name} lacks metadata samples: {absent[:5]}')
				continue
			fields = line.rstrip('\r\n').split('\t')
			chrom, pos_text, ref = fields[0], fields[1], fields[3]
			alts = fields[4].split(',')
			if not is_sequence_allele(ref) or not any(is_sequence_allele(alt) and alt != ref for alt in alts):
				continue
			qual = number_or_none(fields[5], float)
			if qual is None or qual < args.min_site_qual:
				continue
			pos = int(pos_text)
			start = pos - 1
			genes = overlapping_genes(intervals, chrom, start, start + len(ref))
			if not genes:
				continue
			format_indices = {key: index for index, key in enumerate(fields[8].split(':'))}
			parsed = [
				parse_passing_genotype(
					fields[column_index[sample]].split(':'), format_indices, args.min_dp, args.min_gq
				)
				for sample in sample_ids
			]

			for allele_index, alt in enumerate(alts, start=1):
				if not is_sequence_allele(alt) or alt == ref:
					continue
				carriers = np.zeros(len(sample_ids), dtype=bool)
				reference = np.zeros(len(sample_ids), dtype=bool)
				tested = np.zeros(len(sample_ids), dtype=bool)
				depths = []
				qualities = []
				for sample_index, call in enumerate(parsed):
					if call is None:
						continue
					genotype, depth, quality = call
					if allele_index in genotype:
						carriers[sample_index] = True
					elif all(allele == 0 for allele in genotype):
						reference[sample_index] = True
					else:
						continue
					tested[sample_index] = True
					if depth is not None:
						depths.append(depth)
					if quality is not None:
						qualities.append(quality)
				n_carriers = int(carriers.sum())
				n_reference = int(reference.sum())
				call_rate = float(tested.mean())
				if call_rate < args.min_call_rate or min(n_carriers, n_reference) < args.min_carriers:
					continue
				row = {
					'CHROM': chrom, 'POS': pos, 'REF': ref, 'ALT': alt,
					'VARIANT_TYPE': variant_type(ref, alt),
					('RE' + 'F_LENGTH'): len(ref), 'ALT_LENGTH': len(alt),
					'INDEL_LENGTH': len(alt) - len(ref),
					'GENES': ';'.join(genes), 'SITE_QUAL': qual,
					'N_SAMPLES': len(sample_ids), 'N_TESTABLE': int(tested.sum()),
					'N_CARRIERS': n_carriers, 'N_HOM_REF': n_reference,
					'N_OTHER_ALT_OR_LOW_QUALITY': len(sample_ids) - int(tested.sum()),
					'CALL_RATE': call_rate,
					'MEDIAN_DP': float(np.median(depths)) if depths else np.nan,
					'MEDIAN_GQ': float(np.median(qualities)) if qualities else np.nan,
				}
				for phenotype in PHENOTYPES:
					tables = allele_tables(carriers, reference, tested, phenotypes[phenotype], cohorts)
					counts = count_table(tables)
					row[f'{phenotype}_NEG_REF'] = counts[0]
					row[f'{phenotype}_POS_REF'] = counts[1]
					row[f'{phenotype}_NEG_ALT'] = counts[2]
					row[f'{phenotype}_POS_ALT'] = counts[3]
					row[f'{phenotype}_N_TESTED'] = sum(counts)
					row[f'{phenotype}_P_CMH'] = cmh_pvalue(tables)
				rows.append(row)
	return rows


def main():
	parser = argparse.ArgumentParser(
		description='Run one exon-overlapping sequence-allele CMH test family across five binary phenotypes.'
	)
	parser.add_argument('--vcf-dir', required=True, type=Path)
	parser.add_argument('--meta', required=True, type=Path)
	parser.add_argument('--protein-coding-exons', required=True, type=Path)
	parser.add_argument('--outdir', required=True, type=Path)
	parser.add_argument('--chromosomes', nargs='+', choices=CHROMOSOMES, default=CHROMOSOMES)
	parser.add_argument('--min-dp', type=int, default=10)
	parser.add_argument('--min-gq', type=float, default=20)
	parser.add_argument('--min-site-qual', type=float, default=20)
	parser.add_argument('--min-call-rate', type=float, default=0.8)
	parser.add_argument('--min-carriers', type=int, default=2)
	args = parser.parse_args()
	if args.outdir.exists():
		parser.error(f'Output directory already exists: {args.outdir}')
	if not 0 <= args.min_call_rate <= 1 or args.min_carriers < 1:
		parser.error('Call rate must be in [0, 1] and minimum carrier count must be positive')
	missing_vcfs = [args.vcf_dir / f'hs1.{chrom}.vcf' for chrom in args.chromosomes if not (args.vcf_dir / f'hs1.{chrom}.vcf').is_file()]
	if missing_vcfs:
		parser.error('Missing chromosome VCFs: ' + ', '.join(path.name for path in missing_vcfs))
	if not args.protein_coding_exons.is_file() or not args.meta.is_file():
		parser.error('Metadata and protein-coding exon BED must exist')

	metadata = pd.read_csv(args.meta, index_col=0)
	missing_columns = [name for name in PHENOTYPES + ['Cohort'] if name not in metadata]
	if missing_columns:
		parser.error('Metadata lacks columns: ' + ', '.join(missing_columns))
	sample_ids = metadata.index.astype(str).tolist()
	cohorts = metadata['Cohort'].fillna('').astype(str).to_numpy()
	phenotypes = {name: metadata[name].fillna('').astype(str).to_numpy() for name in PHENOTYPES}
	intervals = load_coding_exon_intervals(args.protein_coding_exons)
	rows = []
	for chrom in args.chromosomes:
		vcf_path = args.vcf_dir / f'hs1.{chrom}.vcf'
		chrom_rows = scan_vcf(vcf_path, sample_ids, cohorts, phenotypes, intervals, args)
		rows.extend(chrom_rows)
		print(f'{chrom}: retained {len(chrom_rows)} allele records', flush=True)
	if not rows:
		parser.error('No allele records passed the specified filters')

	results = pd.DataFrame(rows)
	p_columns = [f'{phenotype}_P_CMH' for phenotype in PHENOTYPES]
	pvalues = results[p_columns].to_numpy(dtype=float)
	valid = np.isfinite(pvalues)
	qvalues = np.full(pvalues.shape, np.nan)
	if valid.any():
		qvalues[valid] = multipletests(pvalues[valid], method='fdr_bh')[1]
	for index, phenotype in enumerate(PHENOTYPES):
		results[f'{phenotype}_Q_BH_ALL_TESTS'] = qvalues[:, index]

	args.outdir.mkdir(parents=True)
	output_path = args.outdir / 'all_tests.csv.gz'
	results.to_csv(output_path, index=False, compression='gzip')
	valid_test_count = int(valid.sum())
	readme = f'''# Exon-Overlapping Sequence-Allele CMH Scan

## Test family

Each hypothesis is one alternate sequence allele at one CHM13v2.0 VCF record, tested against one binary phenotype. The five phenotypes are SjD Diagnosis, SSA, SSB, ANA, and RF. Only one statistic is run: a dominant allele-carrier Cochran-Mantel-Haenszel test stratified by cohort. Genotypes carrying a different ALT at the same multiallelic record are excluded from that allele's comparison, not counted as reference.

All finite CMH p-values across every retained allele and all five phenotypes ({valid_test_count:,} tests) are corrected together once using Benjamini-Hochberg. The `*_Q_BH_ALL_TESTS` columns contain parts of this same pooled correction; there are no Fisher or missingness discovery tests in this output.

## Variant and call filters

- Protein-coding, immunoglobulin, or T-cell-receptor transcript-exon overlap in `{args.protein_coding_exons.name}`; a deletion/replacement is retained if its VCF REF span overlaps an exon.
- Literal A/C/G/T sequence alleles only. SNPs, MNVs, insertions, and deletions are included with no maximum allele-length cutoff; symbolic alleles and `*` spanning-deletion alleles are excluded.
- Site QUAL >= {args.min_site_qual:g}; sample DP >= {args.min_dp}; sample GQ >= {args.min_gq:g}.
- For each focal ALT allele, minimum testable call rate >= {args.min_call_rate:.0%} and at least {args.min_carriers} allele carriers and {args.min_carriers} homozygous-reference samples. Other-ALT genotypes and failed/missing calls are not testable.

## Interpretation limits

The VCF alleles are tested in their source representation and were not left-normalized because a matching reference FASTA was not available. Repeated-equivalent representations may therefore remain separate hypotheses. Exon overlap is not a coding-sequence or protein-consequence annotation. These RNA-derived calls are exploratory and are not independently confirmed germline genotypes.
'''
	(args.outdir / 'README.md').write_text(readme, encoding='utf-8')
	print(f'Wrote {len(results)} allele records and {valid_test_count} pooled CMH tests to {output_path}')
	for threshold in (0.05, 0.10, 0.25):
		count = int((qvalues < threshold).sum())
		print(f'Pooled BH q < {threshold:g}: {count} allele-phenotype tests')


if __name__ == '__main__':
	main()