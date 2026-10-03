import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


PHENOTYPES = ['SjD Diagnosis', 'SSA', 'SSB', 'ANA', 'RF']
P_VALUE_COLUMNS = [
	(f'{phenotype}_P_FISHER', f'{phenotype}_Q_FISHER')
	for phenotype in PHENOTYPES
] + [
	(f'{phenotype}_P_CMH', f'{phenotype}_Q_CMH')
	for phenotype in PHENOTYPES
] + [
	(f'{phenotype}_P_MISSING_CMH', f'{phenotype}_Q_MISSING_CMH')
	for phenotype in PHENOTYPES
]


def main():
	parser = argparse.ArgumentParser(
		description='Combine per-chromosome VCF association scans and apply genome-wide BH correction.'
	)
	parser.add_argument('--scan-dir', required=True, type=Path)
	parser.add_argument('--out', required=True, type=Path)
	parser.add_argument('--pattern', default='chr*.assoc.csv.gz')
	args = parser.parse_args()
	if args.out.exists():
		parser.error(f'Output already exists; choose a new path: {args.out}')
	files = sorted(args.scan_dir.glob(args.pattern))
	if not files:
		parser.error(f'No scan files matched {args.scan_dir / args.pattern}')

	frames = [pd.read_csv(path) for path in files]
	results = pd.concat(frames, ignore_index=True)
	for p_column, q_column in P_VALUE_COLUMNS:
		if p_column not in results:
			parser.error(f'Scan file is missing required column: {p_column}')
		pvalues = pd.to_numeric(results[p_column], errors='coerce')
		valid = pvalues.notna()
		results[q_column] = np.nan
		if valid.any():
			results.loc[valid, q_column] = multipletests(pvalues[valid], method='fdr_bh')[1]

	cmh_columns = [f'{phenotype}_P_CMH' for phenotype in PHENOTYPES]
	cmh_values = results[cmh_columns].to_numpy(dtype=float)
	valid = np.isfinite(cmh_values)
	if valid.any():
		q_global = multipletests(cmh_values[valid], method='fdr_bh')[1]
		q_matrix = np.full(cmh_values.shape, np.nan)
		q_matrix[valid] = q_global
		for index, phenotype in enumerate(PHENOTYPES):
			results[f'{phenotype}_Q_CMH_ALL_PHENOTYPES'] = q_matrix[:, index]

	args.out.parent.mkdir(parents=True, exist_ok=True)
	results.to_csv(args.out, index=False, compression='gzip')
	print(f'Combined {len(files)} chromosome files and {len(results)} eligible variant sites.')
	print('Q columns are BH-adjusted across the tested genome; Q_CMH_ALL_PHENOTYPES also spans all five phenotype families.')
	print(f'Wrote {args.out}')


if __name__ == '__main__':
	main()