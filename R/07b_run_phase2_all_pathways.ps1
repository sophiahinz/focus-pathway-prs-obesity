# ==============================================================================
# 07b_run_phase2_all_pathways.ps1
# Runs 07_pathway_pipeline.R sequentially for all pathways in
# <FOCUS_BASE_DIR>/Pathway_Selection/Pathway_Liste.csv,
# to completion (run_phase = 2). Requires the VEP annotation file
# <FOCUS_BASE_DIR>/Outputs_Pathway/Pathway_<pathway_name>_VEP.txt for each pathway.
#
# Hinz et al. Pathway-restricted polygenic scores for energy-substrate
# metabolism and obesity (FoCus cohort). See README for inputs and run order.
#
# Usage:
#   .\07b_run_phase2_all_pathways.ps1                          # all pathways
#   .\07b_run_phase2_all_pathways.ps1 -Only Citrate_Cycle      # selected pathways only
#
# Rscript is taken from the PATH unless the environment variable RSCRIPT
# points to Rscript.exe. Logs: <FOCUS_RESULTS_DIR>/Batch_Logs.
# ==============================================================================
param([string[]]$Only)

$BaseDir     = if ($env:FOCUS_BASE_DIR) { $env:FOCUS_BASE_DIR } else { "path/to/FoCus" }
$ResultsDir  = if ($env:FOCUS_RESULTS_DIR) { $env:FOCUS_RESULTS_DIR } else { "path/to/results" }
$RscriptPath = if ($env:RSCRIPT) { $env:RSCRIPT } else { "Rscript" }
$ScriptPath  = Join-Path $PSScriptRoot "07_pathway_pipeline.R"
$PathwayListe = Join-Path $BaseDir "Pathway_Selection/Pathway_Liste.csv"
$LogDir      = Join-Path $ResultsDir "Batch_Logs"

if (-not (Test-Path $LogDir)) {
    New-Item -ItemType Directory -Path $LogDir | Out-Null
}

$pathways = Import-Csv -Path $PathwayListe
if ($Only) { $pathways = $pathways | Where-Object { $Only -contains $_.pathway_name } }

$summary = @()

foreach ($p in $pathways) {

    $name = $p.pathway_name
    $id   = $p.pathway_selection_id
    $logFile = Join-Path $LogDir "$($name)_Phase2.log"

    Write-Host "==> Phase 2: $name ($id)" -ForegroundColor Cyan

    & $RscriptPath $ScriptPath $name $id "2" *> $logFile
    $exitCode = $LASTEXITCODE

    if ($exitCode -eq 0) {
        Write-Host "    OK - Log: $logFile" -ForegroundColor Green
        $status = "OK"
    } else {
        Write-Host "    ERROR (Exit $exitCode) - Log: $logFile" -ForegroundColor Red
        $status = "ERROR (Exit $exitCode)"
    }

    $summary += [PSCustomObject]@{
        pathway_name = $name
        pathway_id   = $id
        status       = $status
        log          = $logFile
    }
}

$summaryFile = Join-Path $LogDir "Phase2_Summary.csv"
$summary | Export-Csv -Path $summaryFile -NoTypeInformation -Encoding UTF8

Write-Host ""
Write-Host "Phase 2 completed for $(@($pathways).Count) pathway(s)." -ForegroundColor Yellow
Write-Host "Summary: $summaryFile"

