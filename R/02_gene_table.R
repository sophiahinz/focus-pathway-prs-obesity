# ==============================================================================
# 02_gene_table.R
# Builds the pathway-gene table Pathway_Gene_Selection.xlsx (one row per gene
# per pathway) from the gene sets of 01_pathway_selection.R.
#
# Hinz et al. Pathway-restricted polygenic scores for energy-substrate
# metabolism and obesity (FoCus cohort). See README for inputs and run order.
# ==============================================================================

#======================================================================================
# Column format:
#   Pathway | Protein | Gene | Chromosome region | KEGG | Reactome |
#   Datensatz-Art | Core Enzyme | Alternativer Gen-Name
# "Protein" and "Chromosome region" are filled by 03_gene_annotation_grch37.R.
#
# "KEGG erweitert" (KEGG extended; two-zone overlap model): for pathways with
# a "KEGG erweitert" variant only the extended, complete set is used (with
# Datensatz-Art = "KEGG erweitert" and the KEGG/Reactome flag according to
# the actual gene origin); pathways without an extended variant are used
# with their original Datensatz-Art (KEGG or Reactome).
#
# The script stops if Pathway_Gene_Selection.xlsx already exists, so that an
# existing file is never overwritten.
#======================================================================================

#======================================================================================
# PACKAGES
#======================================================================================
library(dplyr)
library(stringr)
library(openxlsx)

#======================================================================================
# INPUT/OUTPUT FILES
#======================================================================================
export_dir <- file.path(Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus"), "Pathway_Selection")

supplement_csv <- file.path(export_dir, "Pathway_Selection_Supplement_GeneSets.csv")
excel_file     <- file.path(export_dir, "Pathway_Gene_Selection.xlsx")

# Safety net: do NOT silently overwrite an existing file
if (file.exists(excel_file)) {
  stop(
    "Pathway_Gene_Selection.xlsx already exists at:\n", excel_file,
    "\nPlease first check whether it already contains manual entries (Core Enzyme, ",
    "Alternativer Gen-Name) that would be lost on regeneration. ",
    "If a rebuild is explicitly wanted, deliberately remove this stop ",
    "or rename/back up the existing file first."
  )
}

#======================================================================================
# READ DATA (long format from the pathway-selection pipeline, phase 7)
#======================================================================================
supplement <- read.csv(supplement_csv, stringsAsFactors = FALSE, encoding = "UTF-8")

# read.csv replaces hyphens in column names with dots - harmonize here
names(supplement)[names(supplement) == "Pathway.ID"]     <- "Pathway_ID"
names(supplement)[names(supplement) == "Pathway.Name"]   <- "Pathway_Name"
names(supplement)[names(supplement) == "Gen.Ursprung"]    <- "Gen_Ursprung"
names(supplement)[names(supplement) == "Datensatz.Art"]   <- "Datensatz_Art"

stopifnot(all(c("Pathway_ID", "Pathway_Name", "Gen", "Gen_Ursprung", "Datensatz_Art") %in%
                names(supplement)))

cat("Read:", nrow(supplement), "rows from", basename(supplement_csv), "\n")

#======================================================================================
# STEP 1: DERIVE BASE PATHWAY NAME (display name without species suffix and
# without "(KEGG erweitert)" label, analogous to the existing file)
#======================================================================================
supplement <- supplement %>%
  mutate(
    Pathway_Base = Pathway_Name %>%
      str_remove(fixed(" - Homo sapiens (human)")) %>%
      str_remove(fixed(" (KEGG erweitert)")) %>%
      str_trim()
  )

#======================================================================================
# STEP 2: PER BASE PATHWAY DECIDE WHETHER THE "ERWEITERT" (EXTENDED) VARIANT OR
# THE ORIGINAL ROWS ARE USED (two-zone model, see above)
#======================================================================================
has_extended_variant <- supplement %>%
  filter(Datensatz_Art == "KEGG erweitert") %>%
  distinct(Pathway_Base) %>%
  pull(Pathway_Base)

cat(length(has_extended_variant),
    " pathway(s) with 'KEGG erweitert' variant found:\n", sep = "")
if (length(has_extended_variant) > 0) {
  cat(paste(" -", has_extended_variant), sep = "\n")
}

final_rows <- supplement %>%
  filter(
    (Pathway_Base %in% has_extended_variant & Datensatz_Art == "KEGG erweitert") |
    (!(Pathway_Base %in% has_extended_variant) & Datensatz_Art %in% c("KEGG", "Reactome"))
  )

# Check: per Pathway_Base + gene there must be exactly one row
dupl_check <- final_rows %>%
  count(Pathway_Base, Gen) %>%
  filter(n > 1)

if (nrow(dupl_check) > 0) {
  warning(nrow(dupl_check), " duplicate pathway-gene combinations found - please check:")
  print(dupl_check)
}

#======================================================================================
# STEP 3: DERIVE KEGG/REACTOME FLAG (X) PER GENE FROM GEN_URSPRUNG
#======================================================================================
final_rows <- final_rows %>%
  mutate(
    KEGG     = ifelse(Gen_Ursprung == "KEGG", "X", NA_character_),
    Reactome = ifelse(Gen_Ursprung == "Reactome", "X", NA_character_)
  )

#======================================================================================
# STEP 4: BUILD TARGET TABLE IN THE EXISTING COLUMN FORMAT
#======================================================================================
gene_table <- final_rows %>%
  transmute(
    Pathway = Pathway_Base,
    Protein = NA_character_,              # filled by 03_gene_annotation_grch37.R
    Gene = Gen,
    `Chromosome region` = NA_character_,  # filled by 03_gene_annotation_grch37.R
    KEGG = KEGG,
    Reactome = Reactome,
    `Datensatz-Art` = Datensatz_Art,
    `Core Enzyme` = NA_character_,        # remains manual/editorial
    `Alternativer Gen-Name` = NA_character_  # remains manual/editorial
  ) %>%
  arrange(Pathway, Gene)

cat("\nTarget table:", nrow(gene_table), "rows,",
    length(unique(gene_table$Gene)), "unique genes,",
    length(unique(gene_table$Pathway)), "unique pathways.\n")

#======================================================================================
# STEP 5: WRITE TO EXCEL (yellow = columns filled in later steps or
# editorial columns)
#======================================================================================
wb <- createWorkbook()
addWorksheet(wb, "Pathway_Gene_Selection")

writeData(wb, sheet = 1, x = gene_table)

manual_style <- createStyle(fgFill = "#FFF2CC")  # yellow: to be filled manually
manual_cols  <- which(names(gene_table) %in% c("Core Enzyme", "Alternativer Gen-Name",
                                                 "Protein", "Chromosome region"))

addStyle(
  wb, sheet = 1, style = manual_style,
  rows = 2:(nrow(gene_table) + 1), cols = manual_cols,
  gridExpand = TRUE, stack = TRUE
)

freezePane(wb, sheet = 1, firstRow = TRUE)
setColWidths(wb, sheet = 1, cols = 1:ncol(gene_table), widths = "auto")

saveWorkbook(wb, excel_file, overwrite = FALSE)

cat("\nDone. New file saved:\n", excel_file, "\n", sep = "")
cat("Next step: run 03_gene_annotation_grch37.R to add protein and ",
    "chromosome region.\n", sep = "")
