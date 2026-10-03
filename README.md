# Sjogren RNA-seq Reanalysis - PhD Thesis Part I

This folder is a curated package of the first, unpublished section of a PhD thesis analysis. It gathers the revised variant scan, DESeq2 and MORFR comparisons, results, figures, and the audit of the earlier `res9`/notebook/MORFR figures. It is intended as a clear record of what was analyzed, not as a claim of independent clinical or genetic validation.

This directory is initialized as a standalone Git repository on `main`; no commit or remote has been created.

## Main Findings

- The single pooled sequence-allele CMH family tested 567,899 allele-phenotype hypotheses across five phenotypes. Five tests had BH q<0.05: HLA-DQA1/SSA, a 97-bp IGHG3 deletion/RF, METTL7A/ANA, SLC17A5/SSA, and GSR/SjD diagnosis.
- HLA-DQA1 showed the clearest descriptive variant-expression consistency across cohorts. Other q<0.05 carrier-expression patterns were small, reversed across strata, or limited by sparse reference genotypes.
- The RF MORFR top-1,000 gene set had RISC/miRNA-related GO terms, including MIR376C and MIR494. These are hypothesis-generating findings; they do not establish a causal small-RNA mechanism.
- Kendall rho-squared and fold-safe MORFR gene rankings show strong internal top-k correspondence. MORFR feature sampling uses Kendall weights, so this is not independent validation.
- `all_gene.csv` is already ComBat-seq adjusted for cohort with SjD diagnosis protected. In the adjusted expression matrix, partial R2 was 0.0067 for cohort and 0.0157 for diagnosis when controlling for the other factor.

## Figures

PDFs are vector exports; PNGs are 300 dpi previews.

- `figures/publication_20261002/variant_association_tiers.pdf`
- `figures/publication_20261002/variant_expression_by_genotype.pdf`
- `figures/publication_20261002/matched_gene_evidence_matrix.pdf`
- `figures/publication_20261002/go_enrichment_by_evidence.pdf`
- `figures/publication_20261002/morfr_variant_overlap.pdf`
- `figures/publication_20261002/kendall_morfr_rank_concordance.pdf`
- `figures/publication_20261002/gene_DESeq2_volcano_panels.pdf`
- `figures/publication_20261002/de_gene_counts_per_contrast.pdf`
- `figures/publication_20261002/rna_expression_pca_cohort_diagnosis.pdf`
- `figures/publication_20261002/clinical_measurement_qc.pdf`

See `docs/figure_story_audit.md` for the question each legacy figure was intended to answer, the weakness in that older analysis, and the replacement interpretation. `docs/results_methods_summary.md` contains the detailed methods, captions, q tiers, and caveats.

## Package Contents

- `src/`: VCF association, integration, clinical, MORFR, and figure-generation scripts.
- `results/`: selected q<0.25 variant results, aggregate DE gene tables, GO terms, clinical summaries, Kendall/MORFR concordance, and saved MORFR rankings/performance.
- `figures/publication_20261002/`: final figures from the revised analysis.
- `docs/`: results/methods summary and legacy-figure audit.

## Data and Reproduction

Per-sample human data are deliberately not included. In particular, this package excludes subject metadata/clinical measurements, ComBat-seq gene-count matrices, VCF genotype calls, FASTQ/BAM/Salmon quantifications, and `res9.csv`/`res9.gen` sample flags. The full genome-wide scan table is also omitted; `results/variant_associations_q25.csv` contains only the 16 alleles with at least one pooled q<0.25 result.

Reproducing the complete scan or sample-level expression plots requires access to the controlled/local inputs: CHM13v2.0 `hs1.chr*.vcf` files, the matching exon BED (`all_ens.bed`), `meta.csv`, `all_gene.csv`, per-contrast DESeq2 tables, and the GO annotation files. The archived `txDE_all.R` references original Linux filesystem paths that must be adapted to an authorized data environment. Use the external scan and figure scripts from `src/` with explicit input/output paths; do not copy those controlled inputs into this repository.

The analysis scripts were run with Python 3.13. Install the Python dependencies listed in `requirements.txt`. The DESeq2 export stage requires R with tximport, DESeq2, sva, rtracklayer, and GenomicFeatures. The original external RNA-seq and VCF inputs are not part of this package, so aggregate outputs and figures are included as the reproducible record of this run.

## Human-Data Release Check

Before making this folder public, confirm that the study consent, IRB/data-use agreement, and applicable repository policy permit public release of the derived figures and aggregate genetic/clinical summaries. `variant_expression_by_genotype` displays unlabeled sample-level points grouped by genotype and phenotype, and some candidate summary rows contain allele counts. If those displays are not approved for public release, omit or aggregate them further before pushing. No subject identifiers or per-sample tables are included here.

## Scientific Limitations

Genotypes are RNA-derived and lack matched-DNA confirmation. VCF indels were not left-normalized because a matching reference FASTA was unavailable, so equivalent representations may remain separate tests. The genotype-expression plots are descriptive, not eQTL tests. The clinical candidate analyses are post-selection exploratory. The GO terms are hierarchically redundant and separately corrected within gene-set tests. MORFR performance is whole-model outer-holdout squared Pearson correlation; rank overlap and Kendall tau are internal correspondence metrics, not candidate-level predictive validation or classical variance explained.
