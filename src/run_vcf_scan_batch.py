import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import subprocess
import sys


CHROMOSOMES = [
	'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
	'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
	'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrM', 'chrX', 'chrY',
]


def run_chromosome(chrom, vcf_dir, metadata, outdir, scanner, args):
	vcf = vcf_dir / f'hs1.{chrom}.vcf'
	output = outdir / f'{chrom}.assoc.csv.gz'
	command = [
		sys.executable, '-u', str(scanner),
		'--vcf', str(vcf),
		'--meta', str(metadata),
		'--out', str(output),
		'--min-dp', str(args.min_dp),
		'--min-gq', str(args.min_gq),
		'--min-site-qual', str(args.min_site_qual),
		'--min-call-rate', str(args.min_call_rate),
		'--min-carriers', str(args.min_carriers),
		'--progress-every', str(args.progress_every),
	]
	if args.protein_coding_exons is not None:
		command.extend(['--protein-coding-exons', str(args.protein_coding_exons)])
	if args.biallelic_snps:
		command.append('--biallelic-snps')
	print(f'[{chrom}] starting', flush=True)
	process = subprocess.Popen(
		command,
		stdout=subprocess.PIPE,
		stderr=subprocess.STDOUT,
		text=True,
		encoding='utf-8',
		errors='replace',
		bufsize=1,
	)
	assert process.stdout is not None
	for line in process.stdout:
		print(f'[{chrom}] {line.rstrip()}', flush=True)
	return chrom, process.wait()


def main():
	parser = argparse.ArgumentParser(description='Run the VCF scanner concurrently by chromosome.')
	parser.add_argument('--vcf-dir', required=True, type=Path)
	parser.add_argument('--meta', required=True, type=Path)
	parser.add_argument('--outdir', required=True, type=Path)
	parser.add_argument('--chromosomes', nargs='+', choices=CHROMOSOMES, default=CHROMOSOMES)
	parser.add_argument('--workers', type=int, default=2)
	parser.add_argument('--min-dp', type=int, default=0)
	parser.add_argument('--min-gq', type=float, default=0)
	parser.add_argument('--min-site-qual', type=float, default=0)
	parser.add_argument('--protein-coding-exons', type=Path)
	parser.add_argument('--biallelic-snps', action='store_true')
	parser.add_argument('--min-call-rate', type=float, default=0.8)
	parser.add_argument('--min-carriers', type=int, default=2)
	parser.add_argument('--progress-every', type=int, default=250000)
	args = parser.parse_args()
	if args.workers < 1:
		parser.error('--workers must be positive')
	if args.min_dp < 0 or args.min_gq < 0 or args.min_site_qual < 0:
		parser.error('DP, GQ, and site QUAL thresholds must be nonnegative')
	if args.outdir.exists():
		parser.error(f'Output directory already exists; choose a new path: {args.outdir}')
	missing = [args.vcf_dir / f'hs1.{chrom}.vcf' for chrom in args.chromosomes if not (args.vcf_dir / f'hs1.{chrom}.vcf').is_file()]
	if missing:
		parser.error('Missing chromosome VCFs: ' + ', '.join(path.name for path in missing))
	if not args.meta.is_file():
		parser.error(f'Metadata file does not exist: {args.meta}')
	if args.protein_coding_exons is not None and not args.protein_coding_exons.is_file():
		parser.error(f'Exon BED does not exist: {args.protein_coding_exons}')

	args.outdir.mkdir(parents=True)
	scanner = Path(__file__).with_name('scan_vcf_associations.py')
	failed = []
	with ThreadPoolExecutor(max_workers=args.workers) as executor:
		futures = {
			executor.submit(run_chromosome, chrom, args.vcf_dir, args.meta, args.outdir, scanner, args): chrom
			for chrom in args.chromosomes
		}
		for future in as_completed(futures):
			chrom, return_code = future.result()
			if return_code:
				failed.append(chrom)
				print(f'[{chrom}] FAILED with exit code {return_code}', flush=True)
			else:
				print(f'[{chrom}] complete', flush=True)
	if failed:
		raise SystemExit('Failed chromosomes: ' + ', '.join(failed))
	print(f'All {len(args.chromosomes)} requested chromosomes complete in {args.outdir}', flush=True)


if __name__ == '__main__':
	main()