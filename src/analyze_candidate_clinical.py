import argparse
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

from scan_vcf_associations import PHENOTYPES, parse_gt


CLINICAL_TRAITS = [
	'Focus Score', 'Schirmer', 'Tear Film Breakup Time',
	'Rose Bengal/van Bijsterveld', 'Salivary Flow', 'C3', 'C4', 'IgG',
	'Lactate Dehydrogenase', 'Gamma Fraction', 'B2M Fraction',
]
SOURCE_NOTES = {
	'Focus Score': 'Clinical focus score; verify histology/source definitions across cohorts.',
	'Schirmer': 'make_meta.py normalizes each cohort by its own maximum; compare within-cohort ranks only.',
	'Tear Film Breakup Time': 'make_meta.py normalizes each cohort by its own maximum; compare within-cohort ranks only.',
	'Rose Bengal/van Bijsterveld': 'make_meta.py normalizes each cohort by its own maximum; MSG coverage is much higher.',
	'Salivary Flow': 'BRA and MSG values are converted from separate source tables; verify units before raw-scale pooling.',
	'C3': 'make_meta.py converts BRA G/L by 100; verify against MSG lab units.',
	'C4': 'make_meta.py converts BRA G/L by 100; verify against MSG lab units.',
	'IgG': 'Source fields appear to use mg/dL; verify assay harmonization.',
	'Lactate Dehydrogenase': 'make_meta.py divides the BRA source by 2; verify conversion and reference ranges.',
	'Gamma Fraction': 'Separate source tables; distributional comparability still needs confirmation.',
	'B2M Fraction': 'Separate source units/conversions; verify before interpreting raw-scale effects.',
}


def rank_normalize_within_cohort(values, cohorts):
	values = pd.to_numeric(values, errors='coerce')
	ranks = values.groupby(cohorts).rank(method='average')
	counts = values.groupby(cohorts).transform('count')
	probability = ((ranks - 0.5) / counts).clip(1e-6, 1 - 1e-6)
	return pd.Series(norm.ppf(probability), index=values.index)


def load_candidate_sites(path, q_cutoff):
	candidates = pd.read_csv(path)
	q_columns = [f'{phenotype}_Q_CMH_ALL_PHENOTYPES' for phenotype in PHENOTYPES]
	candidates['Min_Q_All_Phenotypes'] = candidates[q_columns].min(axis=1)
	candidates = candidates[candidates['Min_Q_All_Phenotypes'] <= q_cutoff].copy()
	return candidates.drop_duplicates(['#CHROM', 'POS', 'REF', 'ALT'])


def extract_carriers(vcf_dir, metadata, candidates, min_dp, min_gq):
	sample_ids = metadata.index.astype(str).tolist()
	candidates_by_chrom = defaultdict(dict)
	for _, candidate in candidates.iterrows():
		key = (str(candidate['#CHROM']), int(candidate['POS']), str(candidate['REF']), str(candidate['ALT']))
		candidates_by_chrom[key[0]].setdefault(key[1], []).append(key)

	carriers = {}
	for chrom, by_pos in candidates_by_chrom.items():
		path = vcf_dir / f'hs1.{chrom}.vcf'
		if not path.is_file():
			raise FileNotFoundError(f'Missing VCF: {path}')
			
		column_index = None
		found = set()
		with path.open(encoding='utf-8', errors='replace') as vcf:
			for line in vcf:
				if line.startswith('##'):
					continue
				if line.startswith('#CHROM\t'):
					header = line.rstrip('\n').split('\t')
					column_index = {sample: i + 9 for i, sample in enumerate(header[9:])}
					absent = [sample for sample in sample_ids if sample not in column_index]
					if absent:
						raise ValueError(f'{path.name} lacks metadata sample IDs: {absent[:5]}')
					continue
				if line.startswith('#'):
					continue
				coordinate = line.split('\t', 2)
				pos = int(coordinate[1])
				if pos not in by_pos:
					if pos > max(by_pos):
						break
					continue
				fields = line.rstrip('\r\n').split('\t')
				ref = fields[3]
				alts = fields[4].split(',')
				format_indices = {key: i for i, key in enumerate(fields[8].split(':'))}
				for key in by_pos[pos]:
					if key in found or key[2] != ref or key[3] not in alts:
						continue
					found.add(key)
					alt_index = alts.index(key[3]) + 1
					site_carriers = {}
					for sample in sample_ids:
						call = fields[column_index[sample]].split(':')
						called, passes, _, dosage, depth, quality, ploidy = parse_gt(
							call, format_indices, min_dp, min_gq
						)
						site_carriers[sample] = {
							'Carrier': (dosage > 0) if passes else pd.NA,
							'Pass': passes,
							'GT_called': called,
							'DP': depth,
							'GQ': quality,
							'Dosage': dosage if passes else pd.NA,
							'Ploidy': ploidy if passes else pd.NA,
						}
					carriers[key] = site_carriers
	missing = [key for keys in candidates_by_chrom.values() for positions in keys.values() for key in positions if key not in carriers]
	if missing:
		raise ValueError(f'Could not find {len(missing)} candidate records in VCFs; examples: {missing[:5]}')
	return carriers


def make_measurement_audit(metadata):
	rows = []
	for trait in CLINICAL_TRAITS:
		values = pd.to_numeric(metadata[trait], errors='coerce')
		for cohort in sorted(metadata['Cohort'].dropna().unique()):
			group = values[metadata['Cohort'] == cohort].dropna()
			rows.append({
				'Trait': trait,
				'Cohort': cohort,
				'N_nonmissing': int(group.size),
				'N_missing': int((metadata.loc[metadata['Cohort'] == cohort, trait].isna()).sum()),
				'Median_raw': group.median() if not group.empty else np.nan,
				'Q25_raw': group.quantile(0.25) if not group.empty else np.nan,
				'Q75_raw': group.quantile(0.75) if not group.empty else np.nan,
				'Harmonization_note': SOURCE_NOTES[trait],
			})
	return pd.DataFrame(rows)


def main():
	parser = argparse.ArgumentParser(
		description='Audit and test strict candidate SNP carriers against continuous clinical traits.'
	)
	parser.add_argument('--associations', type=Path, default=Path('coding_exon_snp_scan_20261002/all_associations_verified.csv.gz'))
	parser.add_argument('--evidence', type=Path, default=Path('coding_exon_snp_scan_20261002/exploratory_evidence_q25.csv'))
	parser.add_argument('--vcf-dir', type=Path, default=Path(r'C:\Users\thoma\Dropbox\HelixBiowulf\2021-02-03_GATK'))
	parser.add_argument('--meta', type=Path, default=Path('meta.csv'))
	parser.add_argument('--out-prefix', type=Path, default=Path('coding_exon_snp_scan_20261002/clinical_q25'))
	parser.add_argument('--min-dp', type=int, default=10)
	parser.add_argument('--min-gq', type=float, default=20)
	parser.add_argument('--q-cutoff', type=float, default=0.25)
	args = parser.parse_args()

	for path in (args.associations, args.evidence, args.meta):
		if not path.is_file():
			parser.error(f'Input file does not exist: {path}')
	if not 0 < args.q_cutoff <= 1:
		parser.error('--q-cutoff must be in (0, 1]')

	metadata = pd.read_csv(args.meta, index_col=0)
	missing_traits = [trait for trait in CLINICAL_TRAITS + ['Age', 'Cohort'] if trait not in metadata]
	if missing_traits:
		parser.error('Metadata lacks columns: ' + ', '.join(missing_traits))
	associations = pd.read_csv(args.associations)
	evidence = pd.read_csv(args.evidence)
	candidates = load_candidate_sites(args.associations, args.q_cutoff)
	genes = {
		(tuple(row[key] for key in ['#CHROM', 'POS', 'REF', 'ALT'])): ';'.join(
			sorted(set(evidence.loc[
				(evidence['#CHROM'] == row['#CHROM'])
				& (evidence['POS'] == row['POS'])
				& (evidence['REF'] == row['REF'])
				& (evidence['ALT'] == row['ALT']), 'Gene'
			].dropna().astype(str)) - {''})
		)
		for _, row in candidates.iterrows()
	}
	genotypes = extract_carriers(args.vcf_dir, metadata, candidates, args.min_dp, args.min_gq)
	args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
	audit_path = args.out_prefix.with_name(args.out_prefix.name + '_measurement_audit.csv')
	results_path = args.out_prefix.with_name(args.out_prefix.name + '_associations.csv')
	for path in (audit_path, results_path):
		if path.exists():
			parser.error(f'Output already exists; will not overwrite: {path}')
	make_measurement_audit(metadata).to_csv(audit_path, index=False)

	rows = []
	for key, sample_calls in genotypes.items():
		carrier = pd.Series({sample: call['Carrier'] for sample, call in sample_calls.items()}, dtype='boolean')
		candidate = candidates[
			(candidates['#CHROM'] == key[0])
			& (candidates['POS'] == int(key[1]))
			& (candidates['REF'] == key[2])
			& (candidates['ALT'] == key[3])
		].iloc[0]
		for trait in CLINICAL_TRAITS:
			values = pd.to_numeric(metadata[trait], errors='coerce')
			rank_score = rank_normalize_within_cohort(values, metadata['Cohort'])
			frame = pd.DataFrame({
				'RankScore': rank_score,
				'Carrier': carrier,
				'Cohort': metadata['Cohort'],
				'Age': pd.to_numeric(metadata['Age'], errors='coerce'),
			}).dropna()
			n_carrier = int(frame['Carrier'].sum())
			n_reference = int((~frame['Carrier'].astype(bool)).sum())
			base = {
				'#CHROM': key[0], 'POS': key[1], 'REF': key[2], 'ALT': key[3],
				'Gene': genes[key], 'Trait': trait,
				'Candidate_Q_All_Phenotypes': candidate['Min_Q_All_Phenotypes'],
				'N': int(len(frame)), 'N_Carrier': n_carrier, 'N_Reference': n_reference,
				'N_missing_trait_or_genotype_or_age': int(len(metadata) - len(frame)),
				'Rank_score_within_cohort': True,
			}
			if len(frame) < 30 or min(n_carrier, n_reference) < 5 or frame['RankScore'].nunique() < 2:
				rows.append({**base, 'Beta_Carrier': np.nan, 'SE_HC3': np.nan, 'P': np.nan, 'Status': 'insufficient data'})
				continue
			design = pd.DataFrame({
				'Carrier': frame['Carrier'].astype(int),
				'Age': frame['Age'],
			}, index=frame.index)
			design = pd.concat([design, pd.get_dummies(frame['Cohort'], prefix='Cohort', drop_first=True, dtype=float)], axis=1)
			design = sm.add_constant(design.astype(float), has_constant='add')
			fit = sm.OLS(frame['RankScore'].astype(float), design).fit(cov_type='HC3')
			ci = fit.conf_int().loc['Carrier']
			rows.append({
				**base,
				'Beta_Carrier': float(fit.params['Carrier']),
				'SE_HC3': float(fit.bse['Carrier']),
				'CI95_Lower': float(ci.iloc[0]),
				'CI95_Upper': float(ci.iloc[1]),
				'P': float(fit.pvalues['Carrier']),
				'Status': 'tested',
			})

	results = pd.DataFrame(rows)
	valid = results['P'].notna()
	results['Q_BH_candidate_clinical_family'] = np.nan
	if valid.any():
		results.loc[valid, 'Q_BH_candidate_clinical_family'] = multipletests(
			results.loc[valid, 'P'], method='fdr_bh'
		)[1]
	results.to_csv(results_path, index=False)
	print(f'Candidate loci analyzed: {len(genotypes)}')
	print(f'Clinical trait tests attempted: {len(results)}; tested: {int(valid.sum())}')
	print(f'Wrote measurement audit: {audit_path}')
	print(f'Wrote candidate clinical associations: {results_path}')
	print('Clinical p-values are exploratory, rank-normalized within cohort, and BH-adjusted across all candidate-by-trait models.')


if __name__ == '__main__':
	main()