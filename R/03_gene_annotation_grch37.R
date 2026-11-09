# ==============================================================================
# 03_gene_annotation_grch37.R
# Adds the official gene/protein name and the GRCh37/hg19 gene coordinates to
# the pathway-gene table and writes Pathway_Gene_Selection_final.xlsx.
#
# Hinz et al. Pathway-restricted polygenic scores for energy-substrate
# metabolism and obesity (FoCus cohort). See README for inputs and run order.
# ==============================================================================
#
# Input : Pathway_Gene_Selection.xlsx (02_gene_table.R)
# Output: Pathway_Gene_Selection_final.xlsx
#           column "Protein"           -> official gene/protein name (GENENAME)
#           column "Chromosome region" -> gene coordinates, GRCh37/hg19
#         Pathway_Gene_Selection_Unmatched_Genes.csv (genes without annotation)
#
# Annotation sources (offline, version-fixed Bioconductor packages):
#   1) org.Hs.eg.db       -> GENENAME and Ensembl gene ID per HGNC symbol
#   2) EnsDb.Hsapiens.v75 -> Ensembl release 75, GRCh37/hg19: gene-level
#                            coordinates (chr:start-end). GRCh37 is consistent
#                            with the FoCus genotype data and the VEP
#                            annotation.
#
# Mapping: HGNC symbol -> Ensembl gene ID (org.Hs.eg.db) -> coordinates
# (EnsDb.Hsapiens.v75, by ID). Genes that cannot be mapped are documented in
# a separate CSV.
# ==============================================================================

library(org.Hs.eg.db)
library(AnnotationDbi)
library(EnsDb.Hsapiens.v75)
library(ensembldb)
library(openxlsx)
library(readxl)
library(dplyr)

export_dir <- file.path(Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus"), "Pathway_Selection")

input_file     <- file.path(export_dir, "Pathway_Gene_Selection.xlsx")
output_file    <- file.path(export_dir, "Pathway_Gene_Selection_final.xlsx")
unmatched_file <- file.path(export_dir, "Pathway_Gene_Selection_Unmatched_Genes.csv")

# Header in row 1, data from row 2
DATA_START_ROW <- 2

#======================================================================================
# READ DATA
#======================================================================================
gene_table <- read_excel(input_file, sheet = 1, col_types = "text")
gene_table <- as.data.frame(gene_table, stringsAsFactors = FALSE)

cat("Read:", nrow(gene_table), "rows,",
    length(unique(gene_table$Gene)), "unique genes.\n")

unique_genes <- sort(unique(gene_table$Gene))

#======================================================================================
# STEP 1: HGNC SYMBOL -> OFFICIAL NAME + ENSEMBL GENE ID (org.Hs.eg.db)
#======================================================================================
symbol_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique_genes,
  keytype = "SYMBOL",
  columns = c("SYMBOL", "GENENAME", "ENSEMBL")
)

# Multiple hits (one symbol -> several Ensembl IDs) are reduced to the first
# hit and reported
multi_ensembl <- symbol_map %>%
  filter(!is.na(ENSEMBL)) %>%
  group_by(SYMBOL) %>%
  filter(n() > 1) %>%
  ungroup()

if (nrow(multi_ensembl) > 0) {
  cat("Note: ", length(unique(multi_ensembl$SYMBOL)),
      " symbol(s) with several Ensembl IDs - first ID used:\n", sep = "")
  print(multi_ensembl)
}

symbol_map <- symbol_map %>%
  filter(!is.na(ENSEMBL)) %>%
  distinct(SYMBOL, .keep_all = TRUE)

no_symbol_match <- setdiff(unique_genes, symbol_map$SYMBOL)

if (length(no_symbol_match) > 0) {
  cat(length(no_symbol_match),
      " gene(s) without hit in org.Hs.eg.db (no GENENAME/Ensembl ID):\n", sep = "")
  cat(paste(no_symbol_match, collapse = ", "), "\n")
}

#======================================================================================
# STEP 2: ENSEMBL GENE ID -> GENOMIC COORDINATES (EnsDb.Hsapiens.v75, GRCh37/hg19)
#======================================================================================
edb <- EnsDb.Hsapiens.v75

gene_ranges <- genes(
  edb,
  filter = GeneIdFilter(unique(symbol_map$ENSEMBL)),
  columns = c("gene_id", "gene_name", "gene_biotype")
)

coord_df <- data.frame(
  ENSEMBL = gene_ranges$gene_id,
  chromosome_region = paste0(
    "chr", as.character(GenomeInfoDb::seqnames(gene_ranges)), ":",
    GenomicRanges::start(gene_ranges), "-",
    GenomicRanges::end(gene_ranges)
  ),
  stringsAsFactors = FALSE
)

annotation <- symbol_map %>%
  left_join(coord_df, by = "ENSEMBL")

no_coord_match <- annotation %>%
  filter(is.na(chromosome_region)) %>%
  pull(SYMBOL)

if (length(no_coord_match) > 0) {
  cat(length(no_coord_match),
      " gene(s) with Ensembl ID but without coordinates in EnsDb.Hsapiens.v75:\n", sep = "")
  cat(paste(no_coord_match, collapse = ", "), "\n")
}

#======================================================================================
# DOCUMENT GENES THAT COULD NOT BE ANNOTATED
#======================================================================================
unmatched_genes <- data.frame(
  Gene = c(no_symbol_match, no_coord_match),
  Reason = c(
    rep("No hit in org.Hs.eg.db (SYMBOL -> GENENAME/ENSEMBL)", length(no_symbol_match)),
    rep("Ensembl ID present, but no gene model in EnsDb.Hsapiens.v75", length(no_coord_match))
  ),
  stringsAsFactors = FALSE
)

if (nrow(unmatched_genes) > 0) {
  write.csv2(unmatched_genes, unmatched_file, row.names = FALSE)
  cat("\nGenes without automatic annotation documented in:\n",
      unmatched_file, "\n", sep = "")
}

#======================================================================================
# STEP 3: TRANSFER ANNOTATION INTO THE TABLE
#======================================================================================
annotation_lookup <- annotation %>%
  select(SYMBOL, GENENAME, chromosome_region) %>%
  distinct(SYMBOL, .keep_all = TRUE)

gene_table <- gene_table %>%
  left_join(annotation_lookup, by = c("Gene" = "SYMBOL")) %>%
  mutate(
    Protein = ifelse(
      !is.na(GENENAME) & (is.na(Protein) | Protein == ""),
      GENENAME,
      Protein
    ),
    `Chromosome region` = ifelse(
      !is.na(chromosome_region) & (is.na(`Chromosome region`) | `Chromosome region` == ""),
      chromosome_region,
      `Chromosome region`
    )
  ) %>%
  select(-GENENAME, -chromosome_region)

cat("\nColumn Protein filled:",
    sum(!is.na(gene_table$Protein) & gene_table$Protein != ""), "of", nrow(gene_table), "rows.\n")
cat("Column Chromosome region filled:",
    sum(!is.na(gene_table[["Chromosome region"]]) & gene_table[["Chromosome region"]] != ""),
    "of", nrow(gene_table), "rows.\n")

#======================================================================================
# STEP 4: WRITE Pathway_Gene_Selection_final.xlsx (formatting of the input
# file is kept; only "Protein" and "Chromosome region" are written)
#======================================================================================
wb <- loadWorkbook(input_file)

protein_col <- which(names(gene_table) == "Protein")
chrom_col   <- which(names(gene_table) == "Chromosome region")

writeData(wb, sheet = 1, x = gene_table$Protein,
          startCol = protein_col, startRow = DATA_START_ROW, colNames = FALSE)
writeData(wb, sheet = 1, x = gene_table[["Chromosome region"]],
          startCol = chrom_col, startRow = DATA_START_ROW, colNames = FALSE)

saveWorkbook(wb, output_file, overwrite = TRUE)

cat("\nDone. File saved:\n", output_file, "\n", sep = "")
cat("org.Hs.eg.db version:      ", as.character(packageVersion("org.Hs.eg.db")), "\n")
cat("EnsDb.Hsapiens.v75 version:", as.character(packageVersion("EnsDb.Hsapiens.v75")),
    "(Ensembl release 75, GRCh37/hg19)\n")
