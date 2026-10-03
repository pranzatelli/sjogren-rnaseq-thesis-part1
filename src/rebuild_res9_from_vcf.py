import argparse
import csv
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


PHENOTYPES = ['SjD Diagnosis', 'SSA', 'SSB', 'ANA', 'RF']


def optional_number(value, cast):
	if value in ('', '.', 'NA'):
		return None
	try:
		return cast(value)
	except ValueError:
		return None


def parse_call(format_keys, sample_value, alt_index, min_dp, min_gq):
	values = sample_value.split(':')
	call = dict(zip(format_keys, values))
	genotype = call.get('GT', '.')
	alleles = genotype.replace('|', '/').split('/')
	called = bool(alleles) and all(allele.isdigit() for allele in alleles)
	dosage = sum(int(allele) == alt_index for allele in alleles) if called else None
	depth = optional_number(call.get('DP', ''), int)
	quality = optional_number(call.get('GQ', ''), float)
	allele_depths = call.get('AD', '.')
	depth_values = allele_depths if allele_depths not in ('', '.') else ''
	passes = called and (min_dp == 0 or (depth is not None and depth >= min_dp))
	passes = passes and (min_gq == 0 or (quality is not None and quality >= min_gq))
	return genotype, depth, quality, depth_values, called, passes, dosage


def parse_legacy_flag(value):
	if pd.isna(value):
		return None
	if isinstance(value, bool):
		return value
	value = str(value).strip().lower()
	if value in ('true', '1'):
		return True
	if value in ('false', '0'):
		return False
	return None


def main():
	parser = argparse.ArgumentParser(
		description='Audit res9 candidate genotypes and rebuild candidate-only carrier tests from VCFs.'
	)
	parser.add_argument('--vcf-dir', required=True, type=Path)
	parser.add_argument('--meta', required=True, type=Path)
	parser.add_argument('--candidates', default=Path('res9.csv'), type=Path)
	parser.add_argument('--outdir', required=True, type=Path)
	parser.add_argument('--min-dp', default=0, type=int)
	parser.add_argument('--min-gq', default=0, type=float)
	args = parser.parse_args()

	candidates = pd.read_csv(args.candidates, index_col=0)
	metadata = pd.read_csv(args.meta, index_col=0)
	missing = [column for column in ('#CHROM', 'POS', 'REF', 'ALT', 'Test', 'Q') if column not in candidates]
	if missing:
		parser.error('Candidate table lacks required columns: ' + ', '.join(missing))
	missing = [column for column in PHENOTYPES if column not in metadata]
	if missing:
		parser.error('Metadata lacks required phenotypes: ' + ', '.join(missing))
	if args.min_dp < 0 or args.min_gq < 0:
		parser.error('--min-dp and --min-gq must be nonnegative')

	sample_ids = metadata.index.astype(str).tolist()
	missing_flags = [sample for sample in sample_ids if sample not in candidates.columns]
	if missing_flags:
		parser.error(f'{len(missing_flags)} metadata samples are absent from candidate columns')

	legacy = {}
	legacy_q = {}
	chrom_keys = {}
	for _, row in candidates.iterrows():
		chrom = str(row['#CHROM'])
		pos = int(row['POS'])
		ref = str(row['REF'])
		alt = str(row['ALT'])
		key = (chrom, pos, ref, alt)
		legacy[key] = {sample: parse_legacy_flag(row[sample]) for sample in sample_ids}
		legacy_q[(key, str(row['Test']))] = row['Q']
		chrom_keys.setdefault(chrom, {})[(pos, ref, alt)] = key

	genotype_rows = []
	association_rows = []
	site_rows = []
	found_total = 0
	for chrom, target_positions in chrom_keys.items():
		vcf_path = args.vcf_dir / f'hs1.{chrom}.vcf'
		if not vcf_path.is_file():
			raise FileNotFoundError(f'Missing VCF for candidate chromosome: {vcf_path}')
		wanted = {}
		for pos, ref, alt in target_positions:
			wanted.setdefault(pos, []).append((ref, alt, target_positions[(pos, ref, alt)]))
		found = set()
		max_pos = max(wanted)
		with vcf_path.open(encoding='utf-8', errors='replace') as vcf:
			for line in vcf:
				if line.startswith('##'):
					continue
				if line.startswith('#CHROM\t'):
					header = line.rstrip('\n').split('\t')
					vcf_samples = header[9:]
					sample_to_column = {sample: i + 9 for i, sample in enumerate(vcf_samples)}
					absent = [sample for sample in sample_ids if sample not in sample_to_column]
					if absent:
						raise ValueError(f'{vcf_path.name} lacks metadata samples: {absent[:5]}')
					continue
				if line.startswith('#'):
					continue
				coordinate_fields = line.split('\t', 2)
				pos = int(coordinate_fields[1])
				if pos not in wanted:
					if pos > max_pos:
						break
					continue
				fields = line.rstrip('\r\n').split('\t')
				ref = fields[3]
				alts = fields[4].split(',')
				for target_ref, target_alt, key in wanted[pos]:
					if target_ref != ref or target_alt not in alts:
						continue
					if key in found:
						raise ValueError(f'Duplicate matching VCF record for {key}')
					found.add(key)
					alt_index = alts.index(target_alt) + 1
					format_keys = fields[8].split(':')
					calls = {}
					for sample in sample_ids:
						call = parse_call(
							format_keys,
							fields[sample_to_column[sample]],
							alt_index,
							args.min_dp,
							args.min_gq,
						)
						calls[sample] = call
						genotype, depth, quality, ad, called, passes, dosage = call
						genotype_rows.append({
							'#CHROM': chrom, 'POS': pos, 'REF': ref, 'ALT': target_alt,
							'Sample': sample, 'GT': genotype, 'DP': depth, 'GQ': quality,
							'AD': ad, 'Called': called, 'PassesFilters': passes,
							'ALT_Dosage': dosage,
							'ALT_Carrier': (dosage > 0) if passes else None,
							'Res9_Flag': legacy[key][sample],
						})

					called_flags = [
						(legacy[key][sample], calls[sample][6] > 0)
						for sample in sample_ids
						if calls[sample][4] and legacy[key][sample] is not None
					]
					passed = [sample for sample in sample_ids if calls[sample][5]]
					site_rows.append({
						'#CHROM': chrom, 'POS': pos, 'REF': ref, 'ALT': target_alt,
						'VCF_ALT_Index': alt_index,
						'N_Metadata_Samples': len(sample_ids),
						'N_Called': sum(calls[sample][4] for sample in sample_ids),
						'N_Passing': len(passed),
						'N_Res9_Flag_True': sum(legacy[key][sample] is True for sample in sample_ids),
						'N_Compared': len(called_flags),
						'N_Flag_Mismatches_vs_ALT': sum(flag != carrier for flag, carrier in called_flags),
					})

					for phenotype in PHENOTYPES:
						counts = {'Negative_REF': 0, 'Positive_REF': 0, 'Negative_ALT': 0, 'Positive_ALT': 0}
						for sample in passed:
							group = metadata.at[sample, phenotype]
							if group not in ('Negative', 'Positive'):
								continue
							carrier = calls[sample][6] > 0
							counts[f'{group}_ALT' if carrier else f'{group}_REF'] += 1
						table = [
							[counts['Negative_REF'], counts['Positive_REF']],
							[counts['Negative_ALT'], counts['Positive_ALT']],
						]
						pvalue = fisher_exact(table)[1] if min(map(sum, table)) > 0 else 1.0
						association_rows.append({
							'#CHROM': chrom, 'POS': pos, 'REF': ref, 'ALT': target_alt,
							'Test': phenotype, 'P': pvalue,
							'Legacy_Q': legacy_q.get((key, phenotype), pd.NA),
							**counts,
							'N_Analyzed': sum(map(sum, table)),
						})

			found_total += len(found)
		for pos, ref, alt in target_positions:
			key = target_positions[(pos, ref, alt)]
			if key not in found:
				site_rows.append({
					'#CHROM': chrom, 'POS': pos, 'REF': ref, 'ALT': alt,
					'VCF_ALT_Index': pd.NA, 'N_Metadata_Samples': len(sample_ids),
					'N_Called': 0, 'N_Passing': 0,
					'N_Res9_Flag_True': sum(value is True for value in legacy[key].values()),
					'N_Compared': 0, 'N_Flag_Mismatches_vs_ALT': pd.NA,
				})

	if association_rows:
		pvalues = [row['P'] for row in association_rows]
		qvalues = multipletests(pvalues, method='fdr_bh')[1]
		for row, qvalue in zip(association_rows, qvalues):
			row['Q_within_res9_candidates'] = qvalue

	args.outdir.mkdir(parents=True, exist_ok=True)
	outputs = {
		'res9_vcf_genotypes.csv': genotype_rows,
		'res9_vcf_site_audit.csv': site_rows,
		'res9_vcf_associations.csv': association_rows,
	}
	for filename, rows in outputs.items():
		path = args.outdir / filename
		if path.exists():
			raise FileExistsError(f'Will not overwrite existing output: {path}')
		pd.DataFrame(rows).to_csv(path, index=False)
	print(f'Candidate loci: {len(legacy)}; matched VCF records: {found_total}')
	print(f'Wrote outputs to {args.outdir.resolve()}')
	print('Q_within_res9_candidates corrects only the existing candidate-locus by phenotype tests; it is not genome-wide discovery FDR.')


if __name__ == '__main__':
	main()