# Pathway-restricted polygenic scores for energy-substrate metabolism and obesity – analysis code

Analysis code for:

> Hinz S, Geisler C, Rohmann N, Hartmann K, Türk K, Franke A, Hollstein T, Laudes M.
> *Pathway-restricted polygenic scores for energy-substrate metabolism and obesity: a systematic screen in the FoCus cohort.*
> International Journal of Obesity (submitted).

Archived version: [![DOI](https://zenodo.org/badge/1446046879.svg)](https://doi.org/10.5281/zenodo.23960201)

The repository contains the R, PowerShell and Python scripts that produced all results, tables and figures of the article, and the 54 pathway-level gene sets. **No individual-level data are included.** The individual-level FoCus data are not publicly available owing to data-protection regulations and the terms of participant consent; they can be made available by the FoCus study management (Institute of Diabetes and Clinical Metabolic Research, Kiel) upon reasonable request and subject to an approved data-use agreement.

## Repository structure

```
R/            analysis scripts (R, PowerShell), numbered in run order
python/       table and figure assembly from the R outputs; Figure 1
gene_sets/    pathway gene sets (GMT and long CSV) and the list of the 54 pathways
```

## Run order

| Step | Script | Purpose | Article |
|---|---|---|---|
| 01 | `01_pathway_selection.R` (+ `01_functions_pathway_selection.R`) | Pathway selection: KEGG base list (BRITE), KEGG (Release 119.0) and Reactome (release 97) gene sets, harmonization to HGNC symbols, two-zone Jaccard overlap model, final gene sets and flow counts | Figure 1; Tables S1, S2b |
| 00 | `00_run_preparation.ps1` | Runs steps 02–06 in order | |
| 02 | `02_gene_table.R` | Pathway-gene table from the gene sets | |
| 03 | `03_gene_annotation_grch37.R` | Gene/protein names and GRCh37 gene coordinates (org.Hs.eg.db, EnsDb.Hsapiens.v75) | |
| 04 | `04_pathway_range_files.R` | PLINK range file per pathway | |
| 05 | `05_plink_extraction_qc.ps1` | Per-pathway SNP extraction and QC per recruitment arm (PLINK 1.9) | Suppl. Methods 1 |
| 06 | `06_global_PCA.R` | Genome-wide harmonization diagnostics, QC, LD pruning, KING kinship, PCA (gPC1–gPC3) | Methods; Suppl. Methods 1 |
| 07 | `07_pathway_pipeline.R` (+ `07a`/`07b` batch runners) | Per pathway: merge of the recruitment arms, QC, HWE, LD pruning, train/test split, univariate SNP models, SNP-sum-score, multivariable models and sensitivity variants (unpruned, beta-, ridge-, elastic-net-weighted, functional, LoF, KEGG-only); SKAT (additional analysis, not reported in the article); descriptive statistics | Table 1; Table 2; Tables S3, S7, S8 |
| 08 | `08_testset_replication.R` | Independent test-set replication | Table 2 |
| 09 | `09_aggregation_robust_candidates.R` | FDR correction, robust-candidate criteria, selection flow | Figure 2; Tables S2a, S3 |
| 10 | `10_pathway_gene_overlap_jaccard.R` | Gene overlap between pathways (Jaccard index) | Results |
| 11 | `11_competitive_permutation_test.R` | Competitive test vs. size-matched random SNP sets | Table S11 |
| 12 | `12_testset_or_deltaR2.R` | Test-set OR per SD and ΔR² | Table 2; Figure 3; Results |
| 13a–c | `13a_…prepare.R`, `13b_…score.R`, `13c_…permutation_test.R` | Externally oriented scores (GIANT BMI GWAS) and their permutation test | Results |
| 14a–c | `14a_…cell_counts.R`, `14b_…models.R`, `14c_…maf_sensitivity.R` | SNP-level functional-concentration analysis incl. MAF adjustment | Tables S9, S10 |
| 15a–c | `15a_…age_gam.R`, `15b_…sports_energy.R`, `15c_suppl_figure_S1.R` | Functional form of continuous covariates | Table S4; Figure S1 |
| – | `python/figure1_pathway_selection.py` | Draws Figure 1 | Figure 1 |
| – | `python/build_tables_figures.py` | Assembles Table 2, Tables S2a, S3, S4, S7–S11 and Figures 2–3 from the R outputs | |

Step 01 combines scripted steps with manual review steps (gap search, inclusion/exclusion criteria, final decision per candidate), which were recorded in a selection spreadsheet that is not part of this repository. Its outputs are provided in `gene_sets/`; `00_run_preparation.ps1` copies them to `<FOCUS_BASE_DIR>/Pathway_Selection/` if they are not yet present there, so the analysis can be reproduced from step 02 onwards.

Step 07 runs in two phases per pathway:

1. `07a_run_phase1_all_pathways.ps1` runs each pathway up to LD pruning and writes the pruned SNP set as PLINK files (`Working_Plink_Pathway_<pathway>_after_R_QC_LDpruned.bed/.bim/.fam` in `<FOCUS_RESULTS_DIR>/Working_Plink/plink_pathway_outputs/`).
2. The PLINK files of each pathway are converted to VCF with PLINK 1.9:
   ```
   plink --bfile Working_Plink_Pathway_<pathway>_after_R_QC_LDpruned --recode vcf --out Working_Plink_Pathway_<pathway>_after_R_QC_LDpruned_VCF
   ```
3. Each VCF is annotated with the Ensembl VEP web interface (VEP v116.0, GRCh37; consequence terms, IMPACT, SIFT and PolyPhen as defined in Supplementary Methods 2). The tab-delimited VEP output is downloaded and saved as `<FOCUS_BASE_DIR>/Outputs_Pathway/Pathway_<pathway>_VEP.txt`. Steps 2 and 3 are manual; they are not scripted.
4. `07b_run_phase2_all_pathways.ps1` completes the analysis of each pathway.

The global train/test split (stratified 70/30 by obesity status, seed 123) is drawn once per pathway in step 07 and stored in `<FOCUS_BASE_DIR>/Outputs_Pathway/`; if a split file already exists there, it is reused.

## Configuration

All scripts read their folders from environment variables:

| Variable | Content |
|---|---|
| `FOCUS_BASE_DIR` | Project root with the input data (see below) |
| `FOCUS_PHENO_DIR` | Folder that directly contains the three FoCus phenotype exports (baseline phenotypes, physical activity, nutrition), with the file names expected by the scripts: `BL_pheno_names_clean_CSV_UTF8_neu.csv`, `FO_PA_alldata_CSV_UTF8_neu.csv`, `FO_nutri_mima_CSV_UTF8_neu.csv` |
| `FOCUS_RESULTS_DIR` | Root folder for the results of steps 07–15 (`OutPW/`, `OutExt/`, `Outputs_FunctionalForm/`, `Working_Plink/`, `Batch_Logs/`) |
| `RSCRIPT`, `PLINK` | Optional: paths to `Rscript.exe` and PLINK 1.9 if they are not on the PATH / in `Working_Plink/` |
| `TABLES_OUT_DIR` | Optional: output folder of the Python scripts (default `./tables_figures_output`) |

Example: `Sys.setenv(FOCUS_BASE_DIR = "/path/to/FoCus")` in R, or `$env:FOCUS_BASE_DIR = "C:/path/to/FoCus"` in PowerShell.

Expected sub-folders of `FOCUS_BASE_DIR`:

| Folder | Content |
|---|---|
| `Genetik_Daten_FoCus/` (genotype data) | Genome-wide GSA genotypes per recruitment arm (`FO_GSA_Ambu`, `FO_GSA_MD`; PLINK binary) |
| `Working_Plink/` | `plink.exe` (1.9), `plink2.exe`, `FO_GSA_Ambu`/`FO_GSA_MD`, range files (04) and `plink_pathway_outputs/` (05) |
| `Pathway_Selection/` | Gene sets and pathway list (01 or `gene_sets/`), pathway-gene tables (02, 03) |
| `Global_PCA/` | Genome-wide PCs and pruned genome-wide dataset (06); `high_LD_regions_hg19.txt` |
| `Harmonization_Diagnostics/` | Diagnostic tables of 06 |
| `Outputs_Pathway/` | Global train/test split per pathway (07) and VEP annotation files |
| `Sensitivity_External_SNP_Orientation/` | GIANT GWAS downloads and the external orientation lookup (13a) |

Scripts 07–13 can additionally write run information into an optional methods/results log (`<FOCUS_RESULTS_DIR>/Methods_Results_Log.md`, marker-delimited appendices); if the file does not exist, this step is skipped with a warning.

## Software

- R 4.5.2; packages: snpStats 1.60.0, bigsnpr 1.12.21, glmnet, caret, data.table, dplyr, stringr, readr, KEGGREST, ReactomeContentService4R, org.Hs.eg.db, AnnotationDbi, EnsDb.Hsapiens.v75, ensembldb, GenomicRanges, GenomeInfoDb, openxlsx, readxl, mgcv, lme4, lmerTest, SKAT, ggplot2, pROC, broom, car, pscl, vcd, margins, ggeffects (versions of the core packages: Supplementary Table S6). Bioconductor packages (snpStats, KEGGREST, ReactomeContentService4R, org.Hs.eg.db, AnnotationDbi, EnsDb.Hsapiens.v75, ensembldb, GenomicRanges, GenomeInfoDb) are installed with `BiocManager::install()`.
- PLINK v1.9.0-b.7.9 and PLINK 2 (PC projection, KING kinship)
- Ensembl VEP v116.0 (cache 116_GRCh37, assembly GRCh37.p13)
- Python 3 with matplotlib and openpyxl (table/figure assembly only)
- PowerShell (batch runners; Windows)

## Gene sets

- `gene_sets/Supplement_GeneSets.gmt` – gene sets in GMT format (pathway ID, database name, HGNC symbols)
- `gene_sets/Pathway_Selection_Supplement_GeneSets.csv` – long format (pathway ID, pathway name, gene, gene origin, set type), as read by steps 02, 07 and 10. Set type “KEGG erweitert” denotes KEGG pathways extended by non-overlapping Reactome genes (two-zone model, zone 2).
- `gene_sets/Pathway_Liste.csv` (pathway list) – the 54 pathways analysed (short pathway name, pathway ID, database, macronutrient category)

## Notes

A few column names and category values used as keys in input files are in German (e.g. `Datensatz-Art`, set type; `KEGG erweitert`, KEGG extended).

## License

MIT License (see LICENSE)
