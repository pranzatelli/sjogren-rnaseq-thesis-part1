# RNA-seq Variant and Expression Results

## Main Results

The pooled sequence-allele scan returned five variant-phenotype tests at BH q<0.05: HLA-DQA1/SSA (q=0.00178), an IGHG3 97-bp deletion/RF (q=0.00753), METTL7A/ANA (q=0.00961), SLC17A5/SSA (q=0.03772), and GSR/SjD diagnosis (q=0.03772). PLBD2/ANA had q=0.0700; 11 additional tests fell in the q<0.25 tier. FGD5 was not reproduced among variants passing the coding-exon and genotype filters; no q<0.25 CMH result mapped to FGD5 exon coordinates in the current BED annotation.

Gene-level DESeq2 tables contain 1,536 significant genes for SjD diagnosis, 996 for SSA, 502 for SSB, 1,681 for ANA, and 955 for RF at per-contrast padj<0.05. HLA-DQA1 is upregulated with significant adjusted p-values for SjD diagnosis, SSA, SSB, and ANA; GSR and METTL7A do not show significant gene-level DE for their corresponding variant phenotypes.

In the existing held-out MORFR rankings, HLA-DQA1 ranks 744 for SSA and IGHG3 ranks 795 for RF. Neither rank is frequent across trees (two and one occurrences, respectively). Among variant-phenotype associations at q<0.25, CD74/RF ranks 181; CD74 also ranks 620 for Focus Score, a cross-target result. METTL7A, SLC17A5, GSR, and PLBD2 are not in their matching-target top 1,000. Overall target-level outer-holdout correlation-squared means range from 0.067 for SSB to 0.478 for Age; these scores summarize whole expression models, not individual variant validation.

The corrected GO analysis returned 4,826 term-set results at within-set BH q<0.05: 4,682 from DE gene sets and 144 from MORFR gene sets. No variant-overlap gene set produced a GO term at q<0.05. These are repeated, correlated ontology terms from separate gene-set analyses, not 4,826 independent findings; they do not support a variant-to-pathway claim.

Within each phenotype, significant GO IDs were shared between DE and MORFR as follows: SjD 8/8, SSA 18/18, SSB 12/35, ANA 44/48, RF 34/35. These are mostly broad immune or cellular terms, while DE contributes far larger gene/term sets. In RF MORFR genes, RISC complex is enriched (q=1.53e-4; 12 expressed-background genes), as is mRNA base-pairing translational repressor activity (q=5.94e-8; 11 genes). The expressed RISC set includes MIR376C and MIR494; their RF MORFR ranks are 23 and 21 (26 and 52 tree occurrences), and their RF DESeq2 padj values are 2.53e-7 and 2.98e-13. This is a specific small-RNA hypothesis to test, not evidence of a causal miRNA mechanism.

The Kendall/MORFR comparison supports the project's intended split-aligned ranking hypothesis. Across seven targets, top-k overlaps between all-sample Kendall rho-squared ranks and fold-safe MORFR ranks exceeded the hypergeometric random expectation for all 42 target-by-k tests (BH q<0.05). At k=25, overlap enrichment ranged from 37.6x to 354x; at k=1,000 it ranged from 1.44x to 4.50x. Within the saved MORFR top 1,000, tie-aware Kendall tau-b ranged from 0.158 to 0.282 (BH q<0.05 for all seven targets). This is internal concordance, not independent validation: MORFR feature sampling itself uses training-only Kendall weights, while `var1` uses all samples.

The expression used for PCA and genotype-expression plots is already ComBat-seq adjusted for cohort, with SjD diagnosis supplied as the protected group. On the adjusted log2(CPM+1) matrix, multivariate least-squares partial R² controlling for the other factor was 0.0067 for cohort and 0.0157 for diagnosis (joint R²=0.0228). This supports residual cohort structure being small rather than dominant in this representation; it does not establish that all batch effects are absent. A second pyCombat/ComBat pass is not warranted from these plots alone.

In the q<0.05 variant-expression panels, HLA-DQA1 carriers had higher mean expression in both cohorts and both SSA-status strata (+0.67 to +2.08 log2(CPM+1)). GSR and SLC17A5 differences reversed direction across cohort/status strata, METTL7A differences were small, and IGHG3 had no BRA reference carriers plus only 11 MSG reference samples overall. No inferential genotype-expression test was fitted; these are descriptive comparisons, not eQTL results.

Same-phenotype DE and MORFR gene lists overlapped by 62 genes for SSB to 239 for ANA, about 3.6-4.8x a simple random-set expectation. This supports internal predictive consistency, but both analyses use the same subjects/phenotypes and MORFR feature sampling is Kendall-weighted; it is not independent biological replication.

Within the q<0.05 genotype-expression panels, HLA-DQA1 ALT carriers had higher mean expression in all four cohort-by-SSA-status strata (descriptive differences +0.67 to +2.08 log2(CPM+1)). GSR and SLC17A5 differences changed direction across cohort/status strata; METTL7A differences were small; IGHG3 had no BRA reference carriers and only 11 MSG reference samples overall. No inferential genotype-expression model was fitted, so these are descriptive checks, not eQTL claims.

Same-phenotype DE and MORFR gene lists overlapped by 62 genes for SSB to 239 for ANA, approximately 3.6-4.8x a simple random-set expectation across targets. This supports within-cohort predictive consistency, but both analyses use the same subjects/phenotypes and MORFR feature sampling is Kendall-weighted; it is not independent biological replication.

## Methods

### Sequence-Allele Association

Joint CHM13v2.0 RNA-derived VCFs were scanned for literal A/C/G/T sequence alleles overlapping protein-coding, immunoglobulin, or T-cell-receptor transcript exon blocks. SNVs, MNVs, insertions, and deletions were retained without a maximum allele-length cutoff. Site QUAL>=20 and sample DP>=10/GQ>=20 were required; each focal allele needed >=80% testable calls, at least two carriers, and at least two homozygous-reference calls. At multiallelic records, genotypes carrying another ALT without the focal ALT were excluded rather than counted as reference.

The single association method was a dominant focal-ALT-carrier Cochran-Mantel-Haenszel test stratified by cohort for each of five binary phenotypes (SjD diagnosis, SSA, SSB, ANA, RF). All finite p-values across 113,988 sequence alleles and five phenotypes (567,899 tests) were corrected together once by Benjamini-Hochberg. No Fisher or missingness association family was run in this scan. Transcript-exon overlap is not a protein consequence annotation.

### Differential Expression

Salmon quantifications were imported with tximport. DESeq2 gene-level contrasts used `~ cohort + phenotype`; each contrast's DESeq2 `padj` was interpreted separately. The volcano figure shows per-contrast `padj<0.05`; vertical +/-1 log2 fold-change guides are visual references, not additional selection criteria.

### MORFR and Clinical Follow-up

MORFR feature weights were restricted to training data within each outer split; model scoring used the untouched outer holdout. The saved run used five repeats, 50 trees per repeat, and 1,000 sampled features per tree. Reported performance is squared Pearson correlation (corr^2), not conventional R^2. The tier integration joins variant genes to matching-target and cross-target top-1,000 ranks without rerunning MORFR.

The earlier clinical candidate analysis tested 143 carrier-by-clinical-trait models selected from the strict SNP candidate set, using within-cohort rank normalization, age/cohort adjustment, HC3 standard errors, and BH correction across those tests. Seventeen had q<0.05. These are post-selection exploratory associations, do not include newly added indel candidates, and are not independent validation. Trait units and missingness differ by cohort.

### Genotype-Expression and Cross-Layer Figures

The genotype-expression figure includes only pooled q<0.05 variant-phenotype tests. Each panel shows expression as log2(CPM+1), grouped on x by cohort and focal genotype, and split/color-coded by the matching phenotype's negative/positive status. Box widths and point jitter are separated to show both status groups. Calls pass DP>=10/GQ>=20; other ALT and low-quality calls are excluded. Counts under each x group are negative/positive. These are descriptive comparisons from the same RNA molecules used to call genotype, not eQTL tests or independent evidence.

The cross-layer matrix includes every gene/variant row overlapping a q<0.25 sequence-allele test and reports pooled CMH q for all five phenotypes, including values above 0.25. It shows DESeq2 log2FC only when the matching per-contrast padj<0.05. MORFR ranks are omitted because only three variant-associated genes appeared in matching-target top 1,000; the separate Kendall-MORFR figure is the appropriate rank-comparison view.

ComBat-seq provenance is in `txDE_all.R`: `all_gene.csv` is written from `txi_adj` after `ComBat_seq(counts, batch=cohort, group=disease)`. No Python ComBat step is applied in figure generation.

The rank-concordance figure compares top-k overlap against a hypergeometric expectation over each target's shared non-missing Kendall feature universe. It also reports Kendall tau between the all-sample Kendall rho-squared rank and MORFR rank among genes in the saved MORFR top 1,000. The MORFR split score is a depth-weighted split occurrence (`0.5^node_depth`, averaged across occurrences), not conventional variance explained.

### GO Enrichment

GO over-representation used Goatools and the local human gene2go associations with GO basic 2024-01-17. The background was genes with CPM>=1 in at least 10% of samples, mapped from gene symbols to NCBI Entrez IDs. For each phenotype separately, tested sets were DESeq2 padj<0.05 genes, genes overlapping variant tests with pooled q<0.25, and the top 1,000 fold-safe MORFR genes. Fisher p-values were BH-adjusted within each phenotype/source set. The GO figure selects up to six strongest MORFR terms per phenotype and shows the same term IDs, in the same row order, beside DE q-values. Its color uses `log10(1 + -log10(q))` to compress the q-value range; dot size is the fraction of the source gene set. No variant-candidate gene set had a term at q<0.05.

## Figure Captions

1. **Variant association tiers.** Each point is one sequence-allele/phenotype test passing the pooled BH q<0.25 threshold. The horizontal category labels give the first overlapping gene and `+n` for additional gene overlaps; point shape denotes SNV, MNV, insertion, or deletion. Dashed vertical lines mark q=0.25, 0.10, and 0.05.
2. **MORFR variant overlap.** Four gene-target rows represent three genes: HLA-DQA1/SSA, IGHG3/RF, and CD74/RF plus CD74/Focus Score. Color is each gene's best variant q tier (red q<0.05, blue 0.10<=q<0.25 in this plot); filled circles match the variant phenotype and open squares are cross-target. The right panel shows whole-model outer-holdout corr^2 means +/- SD, not candidate-level validation.
3. **Gene-level differential expression.** Five per-contrast volcano panels use DESeq2 log2 fold-change and `-log10(padj)`. Color marks direction among genes with per-contrast padj<0.05; the dashed horizontal line marks padj=0.05.
4. **Variant expression by genotype.** Five facets show the q<0.05 variant-phenotype tests. Each compares focal-ALT carriers and homozygous-reference samples within BRA and MSG, with narrower boxes/points colored by negative and positive status of the associated phenotype. Counts below x groups list negative/positive sample numbers. These are descriptive same-RNA comparisons; no genotype-expression significance test is implied.
5. **Phenotype-matched evidence matrix.** Rows are q<0.25 variant/gene overlaps; the left panel prints pooled CMH q for all five phenotypes, including q>0.25. The right panel shows all available matching-contrast DESeq2 log2FC values; `*` marks padj<0.05. MORFR ranks are shown separately.
6. **GO enrichment by evidence set.** For each phenotype, up to six strongest MORFR terms are shown in the same row order beside DESeq2 results. Facet headings report full shared-term counts; color is a compressed q transform and dot size is gene-set fraction. No variant-overlap set had a term at within-set q<0.05.
7. **Differential-expression yield.** Counts show genes with per-contrast DESeq2 padj<0.05, alongside the number of genes with a non-missing adjusted p-value for each contrast.
8. **RNA-expression PCA.** PCA of sample-by-gene log2(CPM+1) values from ComBat-seq-adjusted counts, colored separately by cohort and SjD diagnosis; axes report variance explained. Cohort and diagnosis partial R² estimates are provided in the methods text, not on the plot.
9. **Clinical measurement QC.** The left panel shows non-missing fractions by cohort. The right panel shows pairwise rank correlations after within-cohort rank transformation; sample size varies across trait pairs.
10. **Kendall-MORFR rank concordance.** The left heatmap labels observed/expected intersections at six top-k cutoffs; color gives log2 enrichment. The hypergeometric tail tests excess overlap, not model fit. The right panel gives tie-aware Kendall tau-b between all-sample Kendall rho-squared scores and MORFR ranks among saved top-1,000 genes; positive tau means higher Kendall score tends to accompany a smaller/better MORFR rank. MORFR feature sampling is Kendall-enriched, so this is internal correspondence, not independent validation.

## Interpretation Limits

All genotypes are RNA-derived and lack matched-DNA confirmation. VCF indels were not left-normalized because a matching reference FASTA was unavailable, so equivalent representations can remain separate tests. The cohort is modest for genome-wide association and all molecular/clinical layers reuse the same samples. Findings should be treated as exploratory candidates for orthogonal validation.

Publication-ready vector PDFs and 300-dpi PNGs are in `figures/publication_20261002/`. The source-to-story reconstruction and cautions for legacy notebook/MORFR figures are in `figure_story_audit.md`. Legacy PDFs/SVGs were preserved; the notebook has no executed cells and some retained outputs use older exploratory inputs.