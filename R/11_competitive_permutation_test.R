# ==============================================================================
# 11_competitive_permutation_test.R
# Competitive permutation test of the robust candidate pathways against
# size-matched random SNP sets from the genome-wide pool (Supplementary Table S11).
#
# Hinz et al. Pathway-restricted polygenic scores for energy-substrate
# metabolism and obesity (FoCus cohort). See README for inputs and run order.
# ==============================================================================

# ==============================================================================
# B = 1000 permutations per pathway.
#   - Pathway list: read from AllPathways_Results_Tracking_Table.csv
#     (Robust_Candidate == TRUE).
#   - Checkpointing/resume: results are written to the result CSV after each
#     pathway; a new run skips pathways already completed with status "OK"
#     and identical B (FORCE_RERUN <- TRUE forces a complete re-run).
#   - Errors in individual pathways are caught (tryCatch per pathway) and
#     logged with status "ERROR: <message>"; the remaining pathways continue.
#   - Optional: results and session info are written to Appendix E of the
#     methods/results log (markdown), if that file exists.
# ==============================================================================

library(bigsnpr)
library(bigstatsr)
library(dplyr)

# ------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------
output_dir <- file.path(Sys.getenv("FOCUS_RESULTS_DIR", unset = "path/to/results"), "OutPW")
ranges_dir <- file.path(Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus"), "Working_Plink")
genomewide_pcs_file <- file.path(Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus"), "Global_PCA/Global_Genomewide_PCs.rds")
pool_bed_file <- file.path(Sys.getenv("FOCUS_BASE_DIR", unset = "path/to/FoCus"), "Global_PCA/FO_GSA_Merged_QC_Pruned.bed")
methods_log_path <- file.path(Sys.getenv("FOCUS_RESULTS_DIR", unset = "path/to/results"), "Methods_Results_Log.md")
tracking_table_file <- file.path(output_dir, "AllPathways_Results_Tracking_Table.csv")

output_file <- function(...) file.path(output_dir, paste0(...))
ranges_file <- function(pathway_name) {
  file.path(ranges_dir, paste0("pathway_", pathway_name, "_gene_ranges.txt"))
}

checkpoint_file <- output_file("AllPathways_Permutation_Test_Results.csv")

# ------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------
B_FINAL <- 1000
FORCE_RERUN <- FALSE   # TRUE = recompute all pathways, including those already completed

cov_categorical <- c("gender_cat", "smoking_cat", "school_cert_cat", "graduation_cat")
cov_numeric <- c(
  "sample_age", "sports_activity_min_week", "everyday_activity_min_week",
  "energy_kcal", "gPC1", "gPC2", "gPC3"
)

pathway_seed_base <- function(pathway_name) {
  sum(utf8ToInt(pathway_name)) * 100000L
}

# ------------------------------------------------------------------
# Load robust candidate pathways dynamically from the results tracking table
# (single source of truth -- identical to the definition of "robust
# candidates" in Table 2).
# ------------------------------------------------------------------
if (!file.exists(tracking_table_file)) {
  stop("Tracking table not found: ", tracking_table_file)
}

tracking_table <- read.csv2(tracking_table_file, stringsAsFactors = FALSE)

if (!"Robust_Candidate" %in% colnames(tracking_table)) {
  stop("Column 'Robust_Candidate' not found in ", tracking_table_file, ".")
}

is_robust <- tracking_table$Robust_Candidate %in% c("TRUE", "True", "1", "WAHR", TRUE)
PATHWAYS_FINAL <- tracking_table$Pathway[is_robust]

cat("Robust candidate pathways loaded:", length(PATHWAYS_FINAL), "\n")

#==================================================================================
# 1) Load + impute the genome-wide QC-pruned pool once (with cache)
#==================================================================================
pool_rds_file <- sub("\\.bed$", ".rds", pool_bed_file)

if (!file.exists(pool_rds_file)) {
  cat("Creating backing file from", pool_bed_file, "(one-time, takes a while)...\n")
  snp_readBed(pool_bed_file)
}

pool <- snp_attach(pool_rds_file)
map_pool <- pool$map                            # columns checked below
fam_pool <- pool$fam                             # sample ID column resolved below

stopifnot(all(c("chromosome", "physical.pos", "marker.ID") %in% colnames(map_pool)))
map_pool$chromosome <- as.character(map_pool$chromosome)

id_col_fam <- if ("sample.ID" %in% colnames(fam_pool)) "sample.ID" else "IID"
pool_sample_ids <- as.character(fam_pool[[id_col_fam]])

pool_imputed_meta_file <- sub("\\.rds$", "_imputed_meta.rds", pool_rds_file)

if (file.exists(pool_imputed_meta_file)) {
  cat("Loading cached imputed genome-wide pool:", pool_imputed_meta_file, "\n")
  pool_imputed <- readRDS(pool_imputed_meta_file)
  G_pool <- pool_imputed$genotypes
} else {
  cat(
    "Imputing the genome-wide pool once (mean2, ", ncol(pool$genotypes),
    " SNPs) -- may take a few minutes depending on the machine...\n", sep = ""
  )
  G_pool_imputed <- bigsnpr::snp_fastImputeSimple(
    pool$genotypes, method = "mean2", ncores = 1
  )
  pool_imputed <- list(genotypes = G_pool_imputed, map = map_pool, fam = fam_pool)
  saveRDS(pool_imputed, pool_imputed_meta_file)
  G_pool <- G_pool_imputed
  cat("Imputed pool cached at:", pool_imputed_meta_file, "\n")
}

cat(
  "Genome-wide pool ready:", nrow(map_pool), "SNPs x", length(pool_sample_ids),
  "individuals.\n"
)

#==================================================================================
# 2) Helper: reconstruct training covariate frame per pathway

#==================================================================================
build_train_covariate_frame <- function(pathway_name) {

  final_dataset_ld_file <- output_file(
    "Pathway_", pathway_name, "_Combined_QC_LDpruned.rds"
  )
  if (!file.exists(final_dataset_ld_file)) {
    stop("File not found: ", final_dataset_ld_file)
  }
  pathway_data <- readRDS(final_dataset_ld_file)
  pheno_full <- as.data.frame(pathway_data$pheno, stringsAsFactors = FALSE)
  pheno_full$id <- as.character(pheno_full$id)

  split_file <- output_file(
    "Pathway_", pathway_name, "_Global_TrainTest_Split.rds"
  )
  if (!file.exists(split_file)) stop("File not found: ", split_file)
  global_split <- readRDS(split_file)
  train_ids <- as.character(global_split$train_ids)

  df <- pheno_full

  df$gender_num <- as.numeric(as.character(df$gender))
  df$gender_num[df$gender_num %in% c(0, 7777777, 9999999)] <- NA
  df$gender_cat <- factor(df$gender_num, levels = c(1, 2), labels = c("male", "female"))

  df$smoking_habits_3m_num <- as.numeric(as.character(df$smoking_habits_3m))
  df$smoking_habits_3m_num[df$smoking_habits_3m_num %in% c(0, 7777777, 9999999)] <- NA
  df$smoking_cat <- factor(
    case_when(
      df$smoking_habits_3m_num %in% c(1, 2) ~ "never_smoker",
      df$smoking_habits_3m_num == 3 ~ "former_smoker",
      df$smoking_habits_3m_num == 4 ~ "current_smoker",
      TRUE ~ NA_character_
    ),
    levels = c("never_smoker", "former_smoker", "current_smoker")
  )

  df$school_certificate_A_num <- as.numeric(as.character(df$school_certificate_A))
  df$school_certificate_A_num[!df$school_certificate_A_num %in% 1:5] <- NA
  df$school_cert_cat <- factor(
    ifelse(df$school_certificate_A_num %in% c(1, 2, 5), "risk", "reference"),
    levels = c("reference", "risk")
  )

  df$graduation_num <- as.numeric(as.character(df$graduation))
  df$graduation_num[df$graduation_num %in% c(6, 7777777, 9999999)] <- NA
  df$graduation_cat <- factor(
    ifelse(df$graduation_num %in% c(1, 2, 5), "risk", "reference"),
    levels = c("reference", "risk")
  )

  df$sample_age <- as.numeric(as.character(df$sample_age))
  df$sports_activity_min_week <- as.numeric(as.character(df$sports_activity_min_week))
  df$everyday_activity_min_week <- as.numeric(as.character(df$everyday_activity_min_week))
  df$energy_kcal <- as.numeric(gsub(",", ".", as.character(df$energy_kcal)))

  global_gPCs <- readRDS(genomewide_pcs_file)
  global_gPCs$id <- as.character(global_gPCs$id)
  df <- merge(df, global_gPCs, by = "id", all.x = FALSE)

  outlier_rows_all <- c()
  for (v in cov_numeric) {
    x <- df[[v]]
    Q1 <- quantile(x, 0.25, na.rm = TRUE)
    Q3 <- quantile(x, 0.75, na.rm = TRUE)
    IQRv <- IQR(x, na.rm = TRUE)
    outlier_idx <- which(x < Q1 - 1.5 * IQRv | x > Q3 + 1.5 * IQRv)
    outlier_rows_all <- unique(c(outlier_rows_all, outlier_idx))
  }
  if (length(outlier_rows_all) > 0) df <- df[-outlier_rows_all, , drop = FALSE]

  df_train <- df[df$id %in% train_ids, c("id", "obesity", cov_categorical, cov_numeric)]
  df_train <- df_train[complete.cases(df_train), ]

  if (nrow(df_train) < 10) {
    stop("Too few complete training cases for ", pathway_name)
  }

  df_train
}

#==================================================================================
# 3) Helper: observed (real) test statistic from Regression_Results.csv
#==================================================================================
get_observed_statistic <- function(pathway_name) {
  regression_results_file <- output_file(
    "Pathway_", pathway_name,
    "_Multivariate_SNP_SumScore_Regression_Results_01_Main_LDpruned_All_Oriented.csv"
  )
  if (!file.exists(regression_results_file)) {
    stop("File not found: ", regression_results_file)
  }
  reg <- read.csv2(regression_results_file, stringsAsFactors = FALSE)
  reg_ref <- reg[reg$term == "SNP_SumScore", ]
  if (nrow(reg_ref) != 1) stop("SNP_SumScore row not unique in ", regression_results_file)

  train_or_ref <- as.numeric(reg_ref$OR_Score_Per_SD[1])
  list(train_or_ref = train_or_ref, observed_stat = abs(log(train_or_ref)))
}

#==================================================================================
# 4) Helper: permissible sampling pool (genome-wide pool MINUS whole genes of the
#    pathway)
#==================================================================================
get_eligible_pool_indices <- function(pathway_name) {

  rf <- ranges_file(pathway_name)
  if (!file.exists(rf)) stop("Range file not found: ", rf)

  ranges <- read.table(
    rf, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
    col.names = c("CHR", "START", "END", "GENE")
  )
  ranges$CHR <- as.character(ranges$CHR)

  excluded <- rep(FALSE, nrow(map_pool))

  for (chr in unique(ranges$CHR)) {
    chr_ranges <- ranges[ranges$CHR == chr, ]
    chr_idx <- which(map_pool$chromosome == chr)
    if (length(chr_idx) == 0) next
    pos <- map_pool$physical.pos[chr_idx]

    in_any_gene <- rep(FALSE, length(pos))
    for (g in seq_len(nrow(chr_ranges))) {
      in_any_gene <- in_any_gene |
        (pos >= chr_ranges$START[g] & pos <= chr_ranges$END[g])
    }
    excluded[chr_idx[in_any_gene]] <- TRUE
  }

  which(!excluded)
}

#==================================================================================
# 5) Helper: number of SNPs in the actual pathway (count matching)
#==================================================================================
get_n_snps_pathway <- function(pathway_name) {
  f <- output_file("Pathway_", pathway_name, "_Combined_QC_LDpruned.rds")
  if (!file.exists(f)) stop("File not found: ", f)
  d <- readRDS(f)
  ncol(d$genotypes)
}

#==================================================================================
# 6) Core: ONE permutation for a pathway
#==================================================================================
run_one_permutation <- function(
  pathway_name, df_train, eligible_idx, n_snps_draw, row_idx_train_in_pool, seed
) {

  set.seed(seed)
  drawn_cols <- sample(eligible_idx, size = n_snps_draw, replace = FALSE)

  covar_train_mat <- model.matrix(
    reformulate(c(cov_categorical, cov_numeric)), data = df_train
  )[, -1, drop = FALSE]

  univ <- big_univLogReg(
    X = G_pool,
    y01.train = df_train$obesity,
    ind.train = row_idx_train_in_pool,
    ind.col = drawn_cols,
    ncores = 1
  )
  protective <- which(univ$estim < 0)

  geno_train <- G_pool[row_idx_train_in_pool, drawn_cols, drop = FALSE]

  if (anyNA(geno_train)) {
    col_means <- colMeans(geno_train, na.rm = TRUE)
    for (j in seq_len(ncol(geno_train))) {
      na_idx <- is.na(geno_train[, j])
      if (any(na_idx)) geno_train[na_idx, j] <- col_means[j]
    }
  }

  if (length(protective) > 0) {
    geno_train[, protective] <- 2 - geno_train[, protective, drop = FALSE]
  }

  sum_score <- rowSums(geno_train)
  score_sd <- sd(sum_score)
  if (!is.finite(score_sd) || score_sd == 0) return(NA_real_)
  sum_score_sd <- (sum_score - mean(sum_score)) / score_sd

  df_perm <- df_train
  df_perm$SNP_SumScore_SD <- sum_score_sd

  model_formula <- reformulate(
    termlabels = c("SNP_SumScore_SD", cov_categorical, cov_numeric),
    response = "obesity"
  )

  fit <- tryCatch(
    glm(model_formula, data = df_perm, family = binomial),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NA_real_)

  coefs <- summary(fit)$coefficients
  if (!"SNP_SumScore_SD" %in% rownames(coefs)) return(NA_real_)

  abs(coefs["SNP_SumScore_SD", "Estimate"])
}

#==================================================================================
# 7) Full run per pathway: B_FINAL permutations, with error handling
#==================================================================================
run_permtest_for_pathway <- function(pathway_name, B) {

  result_row <- data.frame(
    pathway = pathway_name, Status = NA_character_,
    n_snps_pathway = NA_integer_, n_eligible_pool = NA_integer_,
    B = NA_integer_, n_valid = NA_integer_, observed_stat = NA_real_,
    train_or_ref = NA_real_, null_mean = NA_real_, null_sd = NA_real_,
    p_value = NA_real_, median_time_per_perm_s = NA_real_,
    mean_time_per_perm_s = NA_real_, total_time_s = NA_real_,
    stringsAsFactors = FALSE
  )

  tryCatch({

    df_train <- build_train_covariate_frame(pathway_name)
    row_idx_train_in_pool <- match(df_train$id, pool_sample_ids)
    if (anyNA(row_idx_train_in_pool)) {
      stop(
        "Not all training IDs found in the genome-wide pool (",
        sum(is.na(row_idx_train_in_pool)), " missing)."
      )
    }

    eligible_idx <- get_eligible_pool_indices(pathway_name)
    n_snps_draw <- get_n_snps_pathway(pathway_name)
    observed <- get_observed_statistic(pathway_name)

    if (n_snps_draw > length(eligible_idx)) {
      stop(
        "Pathway has more SNPs (", n_snps_draw, ") than the permissible ",
        "sampling pool (", length(eligible_idx), " SNPs)."
      )
    }

    seed_base <- pathway_seed_base(pathway_name)
    null_stats <- numeric(B)
    timings <- numeric(B)

    for (b in seq_len(B)) {
      t0 <- Sys.time()
      null_stats[b] <- run_one_permutation(
        pathway_name, df_train, eligible_idx, n_snps_draw,
        row_idx_train_in_pool, seed = seed_base + b
      )
      timings[b] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    }

    n_valid <- sum(!is.na(null_stats))
    if (n_valid == 0) stop("All ", B, " permutations failed (NA).")

    p_value <- (1 + sum(null_stats >= observed$observed_stat, na.rm = TRUE)) /
      (1 + n_valid)

    result_row$Status <- "OK"
    result_row$n_snps_pathway <- n_snps_draw
    result_row$n_eligible_pool <- length(eligible_idx)
    result_row$B <- B
    result_row$n_valid <- n_valid
    result_row$observed_stat <- observed$observed_stat
    result_row$train_or_ref <- observed$train_or_ref
    result_row$null_mean <- mean(null_stats, na.rm = TRUE)
    result_row$null_sd <- sd(null_stats, na.rm = TRUE)
    result_row$p_value <- p_value
    result_row$median_time_per_perm_s <- median(timings)
    result_row$mean_time_per_perm_s <- mean(timings)
    result_row$total_time_s <- sum(timings)

  }, error = function(e) {
    result_row$Status <<- paste0("ERROR: ", conditionMessage(e))
  })

  result_row
}

#==================================================================================
# 8) Methods/results log (optional): Appendix E -- upsert + session info
#==================================================================================
log_permutation_to_md <- function(md_path, results_table) {

  if (!file.exists(md_path)) {
    warning("Methods/results log not found: ", md_path, " - not written.")
    return(invisible(NULL))
  }

  log_lines <- readLines(md_path, warn = FALSE, encoding = "UTF-8")
  begin_marker <- which(trimws(log_lines) == "<!-- AUTOLOG-PERM-PILOT:BEGIN -->")
  end_marker <- which(trimws(log_lines) == "<!-- AUTOLOG-PERM-PILOT:END -->")

  if (length(begin_marker) != 1 || length(end_marker) != 1 || end_marker <= begin_marker) {
    warning("Marker for Appendix E not found in ", md_path, " - not written.")
    return(invisible(NULL))
  }

  header_end <- begin_marker + 2
  data_start <- header_end + 1
  data_end <- end_marker - 1
  existing_data_rows <- if (data_start <= data_end) log_lines[data_start:data_end] else character(0)

  safe_cell <- function(x) {
    if (is.null(x) || length(x) == 0 || is.na(x) || x == "") return("")
    gsub("\\|", "/", gsub("[\r\n]+", " ", as.character(x)))
  }
  extract_pathway_from_row <- function(row) {
    cells <- strsplit(row, "\\|")[[1]]
    if (length(cells) < 2) return(NA_character_)
    trimws(cells[2])
  }

  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M")
  new_pathways <- trimws(as.character(results_table$pathway))

  new_rows <- vapply(seq_len(nrow(results_table)), function(i) {
    r <- results_table[i, ]
    if (identical(r$Status, "OK")) {
      sprintf(
        "| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
        safe_cell(r$pathway), timestamp,
        safe_cell(r$n_snps_pathway), safe_cell(r$n_eligible_pool), safe_cell(r$B),
        safe_cell(round(r$observed_stat, 4)), safe_cell(round(r$null_mean, 4)),
        safe_cell(round(r$null_sd, 4)), safe_cell(round(r$p_value, 4)),
        safe_cell(round(r$median_time_per_perm_s, 3)), safe_cell(round(r$total_time_s, 1))
      )
    } else {
      sprintf(
        "| %s | %s | ERROR | | | | | | %s | | |",
        safe_cell(r$pathway), timestamp, safe_cell(r$Status)
      )
    }
  }, character(1))

  if (length(existing_data_rows) > 0) {
    existing_pathways <- vapply(existing_data_rows, extract_pathway_from_row, character(1))
    kept_existing_rows <- existing_data_rows[
      is.na(existing_pathways) | !(existing_pathways %in% new_pathways)
    ]
  } else {
    kept_existing_rows <- character(0)
  }

  final_data_rows <- c(kept_existing_rows, new_rows)
  updated_lines <- c(log_lines[1:header_end], final_data_rows, log_lines[end_marker:length(log_lines)])
  writeLines(updated_lines, md_path, useBytes = TRUE)

  invisible(NULL)
}

log_permutation_sessioninfo_to_md <- function(md_path) {

  if (!file.exists(md_path)) return(invisible(NULL))

  log_lines <- readLines(md_path, warn = FALSE, encoding = "UTF-8")
  begin_marker <- which(trimws(log_lines) == "<!-- AUTOLOG-PERM-PILOT-SESSIONINFO:BEGIN -->")
  end_marker <- which(trimws(log_lines) == "<!-- AUTOLOG-PERM-PILOT-SESSIONINFO:END -->")

  if (length(begin_marker) != 1 || length(end_marker) != 1 || end_marker <= begin_marker) {
    return(invisible(NULL))
  }

  pkg_version <- function(pkg) {
    tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) "not installed")
  }

  info_lines <- c(
    "",
    paste0("- **Last full run:** ", format(Sys.time(), "%Y-%m-%d %H:%M")),
    paste0("- **R version:** ", R.version.string),
    paste0("- **Platform:** ", R.version$platform, " (", Sys.info()[["sysname"]], ")"),
    paste0(
      "- **Package versions:** bigsnpr ", pkg_version("bigsnpr"),
      ", bigstatsr ", pkg_version("bigstatsr"),
      ", dplyr ", pkg_version("dplyr")
    ),
    paste0("- **Script:** 11_competitive_permutation_test.R (B=", B_FINAL, ")"),
    ""
  )

  updated_lines <- c(log_lines[1:begin_marker], info_lines, log_lines[end_marker:length(log_lines)])
  writeLines(updated_lines, md_path, useBytes = TRUE)
  invisible(NULL)
}

#==================================================================================
# 9) Batch run over all candidate pathways, with checkpointing/resume
#==================================================================================
existing_results <- if (file.exists(checkpoint_file) && !FORCE_RERUN) {
  read.csv2(checkpoint_file, stringsAsFactors = FALSE)
} else {
  NULL
}

already_done <- character(0)
if (!is.null(existing_results) && nrow(existing_results) > 0) {
  already_done <- existing_results$pathway[
    existing_results$Status == "OK" & existing_results$B == B_FINAL
  ]
}

if (length(already_done) > 0) {
  cat(
    length(already_done), " of ", length(PATHWAYS_FINAL),
    " pathways already completed with B = ", B_FINAL,
    " (checkpoint) -- skipped:\n  ",
    paste(already_done, collapse = ", "), "\n", sep = ""
  )
}

pathways_to_run <- setdiff(PATHWAYS_FINAL, already_done)

all_results <- if (!is.null(existing_results)) {
  existing_results[!(existing_results$pathway %in% pathways_to_run), ]
} else {
  data.frame()
}

n_total <- length(PATHWAYS_FINAL)
run_start_time <- Sys.time()

for (i in seq_along(pathways_to_run)) {

  pathway_name <- pathways_to_run[i]

  cat("\n==============================\n")
  cat(
    "Pathway ", which(PATHWAYS_FINAL == pathway_name), "/", n_total, ": ",
    pathway_name, "\n", sep = ""
  )
  cat("==============================\n")

  t_pathway_start <- Sys.time()
  row <- run_permtest_for_pathway(pathway_name, B = B_FINAL)
  t_pathway_elapsed <- as.numeric(difftime(Sys.time(), t_pathway_start, units = "secs"))

  all_results <- dplyr::bind_rows(all_results, row)

  # Checkpoint: persist intermediate result immediately (CSV + Appendix E)
  write.csv2(all_results, checkpoint_file, row.names = FALSE)
  log_permutation_to_md(methods_log_path, row)

  if (identical(row$Status, "OK")) {
    cat(
      "  OK - p_value:", round(row$p_value, 4),
      "| observed_stat:", round(row$observed_stat, 3),
      "| null_mean:", round(row$null_mean, 3),
      "| duration:", round(t_pathway_elapsed, 1), "s\n"
    )
  } else {
    cat("  ", row$Status, "\n")
  }

  elapsed_total_min <- as.numeric(difftime(Sys.time(), run_start_time, units = "mins"))
  n_done_this_run <- i
  avg_min_per_pathway <- elapsed_total_min / n_done_this_run
  n_remaining <- length(pathways_to_run) - i
  cat(
    "  Progress: ", n_done_this_run, "/", length(pathways_to_run),
    " (this run) - elapsed: ", round(elapsed_total_min, 1),
    " min - estimated remaining time: ", round(avg_min_per_pathway * n_remaining, 1),
    " min\n", sep = ""
  )
}

log_permutation_sessioninfo_to_md(methods_log_path)

cat("\n====================================================\n")
cat("FULL RUN COMPLETED\n")
print(all_results[, c("pathway", "Status", "n_snps_pathway", "observed_stat",
                       "null_mean", "null_sd", "p_value")])
cat("====================================================\n")
cat(
  "Result file:", checkpoint_file, "\n",
  "Appendix E (methods log) updated:", methods_log_path, "\n",
  sep = ""
)
