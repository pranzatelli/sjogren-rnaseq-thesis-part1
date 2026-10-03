import argparse
import csv
import gzip
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.stats import chi2, fisher_exact


PHENOTYPES = ['SjD Diagnosis', 'SSA', 'SSB', 'ANA', 'RF']


def number_or_none(value, cast):
	if value in ('', '.', 'NA'):
		return None
	try:
		return cast(value)
	except ValueError:
		return None


def parse_gt(sample_fields, format_indices, min_dp, min_gq):
	gt_index = format_indices.get('GT')
	if gt_index is None or gt_index >= len(sample_fields):
		return False, False, False, None, None, None, 0
	gt = sample_fields[gt_index]
	alleles = gt.replace('|', '/').split('/')
	called = bool(alleles) and all(allele.isdigit() for allele in alleles)
	if not called:
		return False, False, False, None, None, None, 0
	dosage = sum(int(allele) > 0 for allele in alleles)
	dp_index = format_indices.get('DP')
	gq_index = format_indices.get('GQ')
	depth = number_or_none(sample_fields[dp_index], int) if dp_index is not None and dp_index < len(sample_fields) else None
	quality = number_or_none(sample_fields[gq_index], float) if gq_index is not None and gq_index < len(sample_fields) else None
	passes = (min_dp == 0 or (depth is not None and depth >= min_dp))
	passes = passes and (min_gq == 0 or (quality is not None and quality >= min_gq))
	return True, passes, dosage > 0, dosage, depth, quality, len(alleles)


def fisher_pvalue(table):
	if any(sum(row) == 0 for row in table) or any(sum(table[row][column] for row in range(2)) == 0 for column in range(2)):
		return 1.0
	return float(fisher_exact(table).pvalue)


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


def make_table(carriers, included, group):
	negative = included & (group == 'Negative')
	positive = included & (group == 'Positive')
	return [
		[int(np.sum(negative & ~carriers)), int(np.sum(positive & ~carriers))],
		[int(np.sum(negative & carriers)), int(np.sum(positive & carriers))],
	]


def load_protein_coding_exons(path):
	intervals = {}
	with path.open(encoding='utf-8') as bed:
		for line in bed:
			fields = line.rstrip('\n').split('\t')
			if len(fields) < 9 or not fields[1].isdigit() or fields[7] != 'protein_coding':
				continue
			start = int(fields[1])
			sizes = [int(value) for value in fields[5].rstrip(',').split(',') if value]
			offsets = [int(value) for value in fields[6].rstrip(',').split(',') if value]
			if len(sizes) != len(offsets):
				raise ValueError(f'Mismatched BED block arrays for {fields[3]}')
			intervals.setdefault(fields[0], []).extend(
				(start + offset, start + offset + size)
				for size, offset in zip(sizes, offsets)
			)

	merged = {}
	for chrom, blocks in intervals.items():
		blocks.sort()
		chrom_blocks = []
		for start, end in blocks:
			if chrom_blocks and start <= chrom_blocks[-1][1]:
				chrom_blocks[-1] = (chrom_blocks[-1][0], max(chrom_blocks[-1][1], end))
			else:
				chrom_blocks.append((start, end))
		merged[chrom] = (
			np.fromiter((block[0] for block in chrom_blocks), dtype=np.int64),
			np.fromiter((block[1] for block in chrom_blocks), dtype=np.int64),
		)
	return merged


def overlaps_intervals(intervals, chrom, zero_based_pos):
	if intervals is None or chrom not in intervals:
		return intervals is None
	starts, ends = intervals[chrom]
	index = np.searchsorted(starts, zero_based_pos, side='right') - 1
	return index >= 0 and zero_based_pos < ends[index]


def main():
	parser = argparse.ArgumentParser(
		description='Stream one joint VCF and write filtered carrier/missingness association summaries.'
	)
	parser.add_argument('--vcf', required=True, type=Path)
	parser.add_argument('--meta', required=True, type=Path)
	parser.add_argument('--out', required=True, type=Path)
	parser.add_argument('--min-dp', type=int, default=0)
	parser.add_argument('--min-gq', type=float, default=0)
	parser.add_argument('--min-site-qual', type=float, default=0)
	parser.add_argument('--protein-coding-exons', type=Path)
	parser.add_argument('--biallelic-snps', action='store_true')
	parser.add_argument('--min-call-rate', type=float, default=0.8, help='Minimum passing-call rate for genotype tests; missingness tests retain lower-call-rate sites.')
	parser.add_argument('--min-carriers', type=int, default=2, help='Minimum passing ALT carriers and REF calls for genotype tests.')
	parser.add_argument('--max-records', type=int, default=0, help='Optional record limit for a smoke test.')
	parser.add_argument('--progress-every', type=int, default=250000, help='Print progress after this many VCF records; 0 disables progress messages.')
	args = parser.parse_args()
	if args.out.exists():
		parser.error(f'Output already exists; choose a new path: {args.out}')
	if args.min_dp < 0 or args.min_gq < 0 or args.min_site_qual < 0 or args.min_carriers < 1:
		parser.error('Depth, quality, and site quality thresholds must be nonnegative; min-carriers must be positive')
	if not 0 <= args.min_call_rate <= 1:
		parser.error('--min-call-rate must be between 0 and 1')

	metadata = pd.read_csv(args.meta, index_col=0)
	missing = [name for name in PHENOTYPES + ['Cohort'] if name not in metadata]
	if missing:
		parser.error('Metadata lacks columns: ' + ', '.join(missing))
	sample_ids = metadata.index.astype(str).tolist()
	cohorts = metadata['Cohort'].fillna('').astype(str).to_numpy()
	groups = {name: metadata[name].fillna('').astype(str).to_numpy() for name in PHENOTYPES}
	if args.protein_coding_exons is not None and not args.protein_coding_exons.is_file():
		parser.error(f'Exon BED does not exist: {args.protein_coding_exons}')
	exon_intervals = load_protein_coding_exons(args.protein_coding_exons) if args.protein_coding_exons else None

	args.out.parent.mkdir(parents=True, exist_ok=True)
	columns = [
		'#CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER', 'N_SAMPLES',
		'N_CALLED', 'N_PASS', 'RAW_CALL_RATE', 'PASS_RATE',
		'N_ALT_CARRIERS', 'N_REF', 'ALT_ALLELE_COUNT', 'ALT_ALLELE_NUMBER',
		'ALT_ALLELE_FREQUENCY', 'MEDIAN_DP', 'MEDIAN_GQ',
	]
	for phenotype in PHENOTYPES:
		columns.extend([
			f'{phenotype}_NEG_REF', f'{phenotype}_POS_REF',
			f'{phenotype}_NEG_ALT', f'{phenotype}_POS_ALT',
			f'{phenotype}_N_TESTED', f'{phenotype}_P_FISHER',
			f'{phenotype}_P_CMH', f'{phenotype}_P_MISSING_CMH',
		])

	read_count = 0
	write_count = 0
	genotype_test_count = 0
	low_call_rate_count = 0
	low_carrier_count = 0
	outside_exon_count = 0
	non_snp_count = 0
	low_site_qual_count = 0
	with args.vcf.open(encoding='utf-8', errors='replace') as vcf, gzip.open(args.out, 'wt', newline='', encoding='utf-8') as output:
		writer = csv.DictWriter(output, fieldnames=columns)
		writer.writeheader()
		for line in vcf:
			if line.startswith('##'):
				continue
			if line.startswith('#CHROM\t'):
				header = line.rstrip('\n').split('\t')
				vcf_samples = header[9:]
				column_index = {sample: index + 9 for index, sample in enumerate(vcf_samples)}
				absent = [sample for sample in sample_ids if sample not in column_index]
				if absent:
					raise ValueError(f'VCF is missing metadata samples: {absent[:5]}')
				continue
			if line.startswith('#'):
				continue
			if args.max_records and read_count >= args.max_records:
				break
			read_count += 1
			if args.progress_every and read_count % args.progress_every == 0:
				print(f'{args.vcf.name}: scanned {read_count} VCF records', file=sys.stderr, flush=True)
			fields = line.rstrip('\r\n').split('\t')
			if args.biallelic_snps and not (
				len(fields[3]) == 1
				and len(fields[4]) == 1
				and fields[3] in 'ACGT'
				and fields[4] in 'ACGT'
			):
				non_snp_count += 1
				continue
			if exon_intervals is not None and not overlaps_intervals(exon_intervals, fields[0], int(fields[1]) - 1):
				outside_exon_count += 1
				continue
			site_quality = number_or_none(fields[5], float)
			if args.min_site_qual and (site_quality is None or site_quality < args.min_site_qual):
				low_site_qual_count += 1
				continue
			format_indices = {key: index for index, key in enumerate(fields[8].split(':'))}
			called = np.zeros(len(sample_ids), dtype=bool)
			passed = np.zeros(len(sample_ids), dtype=bool)
			carriers = np.zeros(len(sample_ids), dtype=bool)
			dosages = np.zeros(len(sample_ids), dtype=int)
			allele_numbers = np.zeros(len(sample_ids), dtype=int)
			depths = []
			qualities = []
			for index, sample in enumerate(sample_ids):
				call = fields[column_index[sample]].split(':')
				is_called, is_passing, is_carrier, dosage, depth, quality, ploidy = parse_gt(
					call, format_indices, args.min_dp, args.min_gq
				)
				called[index] = is_called
				passed[index] = is_passing
				carriers[index] = is_carrier and is_passing
				if is_passing:
					dosages[index] = dosage
					allele_numbers[index] = ploidy
					if depth is not None:
						depths.append(depth)
					if quality is not None:
						qualities.append(quality)
			n_called = int(called.sum())
			n_pass = int(passed.sum())
			n_carriers = int(carriers.sum())
			n_ref = n_pass - n_carriers
			pass_rate = n_pass / len(sample_ids)
			genotype_testable = pass_rate >= args.min_call_rate and min(n_carriers, n_ref) >= args.min_carriers
			if not genotype_testable:
				low_call_rate_count += int(pass_rate < args.min_call_rate)
				low_carrier_count += int(min(n_carriers, n_ref) < args.min_carriers)
				continue

			allele_number = int(allele_numbers.sum())
			allele_count = int(dosages.sum())
			row = {
				'#CHROM': fields[0], 'POS': fields[1], 'REF': fields[3], 'ALT': fields[4],
				'QUAL': fields[5], 'FILTER': fields[6], 'N_SAMPLES': len(sample_ids),
				'N_CALLED': n_called, 'N_PASS': n_pass,
				'RAW_CALL_RATE': n_called / len(sample_ids), 'PASS_RATE': pass_rate,
				'N_ALT_CARRIERS': n_carriers, 'N_REF': n_ref,
				'ALT_ALLELE_COUNT': allele_count,
				'ALT_ALLELE_NUMBER': allele_number,
				'ALT_ALLELE_FREQUENCY': allele_count / allele_number if allele_number else np.nan,
				'MEDIAN_DP': float(np.median(depths)) if depths else np.nan,
				'MEDIAN_GQ': float(np.median(qualities)) if qualities else np.nan,
			}
			for phenotype, group in groups.items():
				known = np.isin(group, ['Negative', 'Positive'])
				eligible = passed & known
				table = make_table(carriers, eligible, group)
				stratified_tables = [make_table(carriers, eligible & (cohorts == cohort), group) for cohort in np.unique(cohorts)]
				missing_tables = []
				for cohort in np.unique(cohorts):
					in_cohort = known & (cohorts == cohort)
					negative = in_cohort & (group == 'Negative')
					positive = in_cohort & (group == 'Positive')
					missing_tables.append([
						[int(np.sum(passed & negative)), int(np.sum(passed & positive))],
						[int(np.sum(~passed & negative)), int(np.sum(~passed & positive))],
					])
				for label, value in zip(
					[
						f'{phenotype}_NEG_REF', f'{phenotype}_POS_REF',
						f'{phenotype}_NEG_ALT', f'{phenotype}_POS_ALT',
					],
					[table[0][0], table[0][1], table[1][0], table[1][1]],
				):
					row[label] = value
				row[f'{phenotype}_N_TESTED'] = sum(map(sum, table))
				row[f'{phenotype}_P_FISHER'] = fisher_pvalue(table) if genotype_testable else np.nan
				row[f'{phenotype}_P_CMH'] = cmh_pvalue(stratified_tables) if genotype_testable else np.nan
				row[f'{phenotype}_P_MISSING_CMH'] = cmh_pvalue(missing_tables) if 0 < n_pass < len(sample_ids) else np.nan
			writer.writerow(row)
			write_count += 1
			genotype_test_count += int(genotype_testable)
	print(f'{args.vcf.name}: scanned {read_count} records; wrote {write_count} genotype/missingness-informative sites, including {genotype_test_count} genotype-tested sites, to {args.out}')
	print(f'Filters: min-DP={args.min_dp}, min-GQ={args.min_gq}, genotype min-call-rate={args.min_call_rate}, genotype min-carriers={args.min_carriers}')
	print(f'Region/allele filters: min-site-QUAL={args.min_site_qual}, protein-coding-exon restriction={args.protein_coding_exons is not None}, biallelic-SNP-only={args.biallelic_snps}')
	print(f'Excluded sites: {low_call_rate_count} below call-rate threshold; {low_carrier_count} below carrier/reference-count threshold')
	print(f'Excluded before genotype parsing: {outside_exon_count} outside protein-coding exons; {non_snp_count} not biallelic A/C/G/T SNPs; {low_site_qual_count} below site-QUAL threshold')


if __name__ == '__main__':
	main()