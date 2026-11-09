# ==============================================================================
# 00_run_preparation.ps1
# Runs the preparation steps 02-06 in order: pathway-gene table, GRCh37 gene
# annotation, PLINK range files, per-pathway PLINK extraction/QC and the
# genome-wide PCA. Afterwards, the pathway pipeline is started with
# 07a_run_phase1_all_pathways.ps1.
#
# Hinz et al. Pathway-restricted polygenic scores for energy-substrate
# metabolism and obesity (FoCus cohort). See README for inputs and run order.
# ==============================================================================
#
# Step 01 (pathway selection) is not run here: it requires the KEGG/Reactome
# web services and the manually maintained selection spreadsheet. Its outputs
# are provided in gene_sets/. If they are not yet present in
# <FOCUS_BASE_DIR>/Pathway_Selection, they are copied from gene_sets/.
#
# Usage:
#   $env:FOCUS_BASE_DIR    = "C:/path/to/FoCus"
#   $env:FOCUS_PHENO_DIR   = "C:/path/to/phenotype_csv"
#   $env:FOCUS_RESULTS_DIR = "C:/path/to/results"
#   .\00_run_preparation.ps1
#
# Rscript is taken from the PATH unless the environment variable RSCRIPT
# points to Rscript.exe. Steps 02/03 are skipped if their output file already
# exists; 04 does not overwrite existing range files. The script stops at the
# first failing step.
# ==============================================================================

$ErrorActionPreference = "Stop"

$BaseDir     = if ($env:FOCUS_BASE_DIR) { $env:FOCUS_BASE_DIR } else { throw "Please set FOCUS_BASE_DIR." }
$RscriptPath = if ($env:RSCRIPT) { $env:RSCRIPT } else { "Rscript" }
$RepoGeneSets = Join-Path (Split-Path $PSScriptRoot -Parent) "gene_sets"
$SelectionDir = Join-Path $BaseDir "Pathway_Selection"

New-Item -ItemType Directory -Force -Path $SelectionDir | Out-Null
foreach ($f in @("Pathway_Liste.csv", "Pathway_Selection_Supplement_GeneSets.csv", "Supplement_GeneSets.gmt")) {
    $target = Join-Path $SelectionDir $f
    if (-not (Test-Path $target)) {
        Copy-Item (Join-Path $RepoGeneSets $f) $target
        Write-Host "Copied gene_sets/$f -> $target"
    }
}

function Invoke-Step($label, [scriptblock]$command) {
    Write-Host ""
    Write-Host "==> $label" -ForegroundColor Cyan
    & $command
    if ($LASTEXITCODE -ne 0) { throw "$label failed (exit code $LASTEXITCODE)." }
}

# 02/03 are skipped if their output already exists (02 never overwrites it)
if (Test-Path (Join-Path $SelectionDir "Pathway_Gene_Selection.xlsx")) {
    Write-Host "02 skipped: Pathway_Gene_Selection.xlsx already exists."
} else {
    Invoke-Step "02 pathway-gene table" { & $RscriptPath (Join-Path $PSScriptRoot "02_gene_table.R") }
}
if (Test-Path (Join-Path $SelectionDir "Pathway_Gene_Selection_final.xlsx")) {
    Write-Host "03 skipped: Pathway_Gene_Selection_final.xlsx already exists."
} else {
    Invoke-Step "03 GRCh37 gene annotation" { & $RscriptPath (Join-Path $PSScriptRoot "03_gene_annotation_grch37.R") }
}
Invoke-Step "04 PLINK range files"       { & $RscriptPath (Join-Path $PSScriptRoot "04_pathway_range_files.R") }
Invoke-Step "05 PLINK extraction + QC"   { & (Join-Path $PSScriptRoot "05_plink_extraction_qc.ps1"); $global:LASTEXITCODE = 0 }
Invoke-Step "06 genome-wide PCA"         { & $RscriptPath (Join-Path $PSScriptRoot "06_global_PCA.R") }

Write-Host ""
Write-Host "Preparation completed. Next step: .\07a_run_phase1_all_pathways.ps1" -ForegroundColor Yellow
