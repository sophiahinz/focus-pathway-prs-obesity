# Pathway-restricted polygenic scores for energy-substrate metabolism and obesity – analysis code

Analysis code for:

> Hinz S, Geisler C, Rohmann N, Hartmann K, Türk K, Franke A, Hollstein T, Laudes M.
> *Pathway-restricted polygenic scores for energy-substrate metabolism and obesity: a systematic screen in the FoCus cohort.*
> International Journal of Obesity (submitted). [DOI to be added]

The repository contains the R scripts (and two Python helper scripts) that produced all results, tables and figures of the article, plus the 54 pathway-level gene sets. **No individual-level data are included.** FoCus genotype and phenotype data are available from the FoCus study management upon reasonable request (see the Data Availability statement of the article).

## Repository structure

```
R/            analysis scripts, numbered in run order
python/       table/figure assembly from the R outputs; Figure 1
gene_sets/    pathway gene sets (GMT and long CSV) and the pathway list used by the pipeline
```

## Run order

| Step | Script | Purpose | Article |
|---|---|---|---|
| 00 | `00_global_PCA.R` | Genome-wide harmonization diagnostics, QC, LD pruning, KING kinship, PCA (gPC1–gPC3) | Methods; Suppl. Methods 1 |
| 01 | `01_pathway_gene_selection.R` | Pathway/gene table from the pathway-selection output | Figure 1; Tables S1, S2b |
| 02 | `02_pathway_gene_overlap_jaccard.R` | Gene overlap between pathways (Jaccard index) | Methods |
| 03 | `03_pathway_pipeline.R` (+ `03a`/`03b` batch runners) | Per-pathway SNP extraction, QC, HWE, LD pruning, train/test split, SNP-sum score, multivariable models and sensitivity variants (unpruned, beta-, ridge-, elastic-net-weighted, functional, LoF) | Table 2; Tables S3, S7, S8 |
| 04 | `04_testset_replication.R` | Independent test-set replication | Table 2 |
| 05 | `05_aggregation_robust_candidates.R` | FDR correction, robust-candidate criteria | Figure 2; Table S2a |
| 06 | `06_competitive_permutation_test.R` | Competitive test vs. size-matched random SNP sets | Table S11 |
| 07 | `07_testset_or_deltaR2.R` | Test-set OR per SD and ΔR² | Results |
| 08a–c | `08a_…prepare.R`, `08b_…score.R`, `08c_…permutation_test.R` | Externally oriented scores (GIANT BMI GWAS) and their permutation test | Results; Figure 3 |
| 09a–c | `09a_…cell_counts.R`, `09b_…models.R`, `09c_…maf_sensitivity.R` | SNP-level functional-concentration analysis incl. MAF adjustment | Tables S9, S10 |
| 10a–c | `10a_…age_gam.R`, `10b_…sports_energy.R`, `10c_suppl_figure_S1.R` | Functional form of continuous covariates | Table S4; Figure S1 |
| – | `python/build_tables_figures.py` | Assembles Table 2, Tables S2a–S11 and Figures 2–3 from the R outputs | |
| – | `python/figure1_pathway_selection.py` | Draws Figure 1 | |

Step 03 runs in two phases per pathway (`03a`: up to VCF export; then external annotation with Ensembl VEP; `03b`: remaining analysis).

## Configuration

All scripts read the project root from the environment variable `FOCUS_BASE_DIR` (e.g. `Sys.setenv(FOCUS_BASE_DIR = "/path/to/FoCus")` in R, or `$env:FOCUS_BASE_DIR = "..."` in PowerShell). Expected sub-folders below this root are those referenced in the scripts (e.g. `Genetik_Daten_FoCus/`, `FoCus Daten/Daten/CSV_repariert/`, `Global_PCA/`, `Pathway_Selection/`, `Working_Plink/`). The folder `Repair/` holds the outputs of the final analysis run, which used the corrected phenotype export (`CSV_repariert`); all results reported in the article stem from this run.

## Software

- R 4.5.2; packages: snpStats 1.60.0, bigsnpr 1.12.21, glmnet, caret, data.table, dplyr, KEGGREST, ReactomeContentService4R, org.Hs.eg.db, AnnotationDbi, openxlsx, readxl, mgcv, lme4, lmerTest, SKAT, ggplot2, pROC, broom, car, pscl (Supplementary Table S6)
- PLINK v1.9.0-b.7.9 and PLINK 2 (PC projection, KING kinship)
- Ensembl VEP v116.0 (cache 116_GRCh37, assembly GRCh37.p13)
- Python 3 with matplotlib (helper scripts only)

## Gene sets

- `gene_sets/Supplement_GeneSets.gmt` – gene sets in GMT format (pathway ID, name, HGNC symbols)
- `gene_sets/Pathway_Selection_Supplement_GeneSets.csv` – long format (pathway ID, name, gene, gene source, set type), as read by `03_pathway_pipeline.R`. Set type “KEGG erweitert” denotes KEGG pathways supplemented with non-overlapping Reactome genes.
- `gene_sets/Pathway_Liste.csv` – the 54 pathways analysed (pathway name, ID, database, macronutrient category)

## Notes

Code comments and console messages are in English. A small number of data labels, column names and status codes (e.g. `Datensatz-Art`, `KEGG erweitert`, `Gen-Ursprung`, internal step labels) are kept in their original German form because the scripts use them as keys when reading and writing intermediate files.

## License

[to be added]
