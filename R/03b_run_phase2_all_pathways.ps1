# ==============================================================================
# REPAIR VERSION (2026-10-09): copy of Run_Phase2_AllPathways.ps1
# Changes: ScriptPath -> _REPAIR.R, LogDir -> Repair/Batch_Logs (marked '# REPAIR').
# Optional: test individual pathways with  .\03b_run_phase2_all_pathways.ps1 -Only Citrate_Cycle
# ==============================================================================
param([string[]]$Only)   # REPAIR

# ==============================================================================
# Run_Phase2_AllPathways.ps1
#
# Runs the R script "03_pathway_pipeline.R"
# sequentially for ALL pathways from "Pathway_Liste.csv" to completion
# (run_phase = 2) - requires that the external
# VCF creation + VEP annotation has already been performed for each pathway and the file
# Pathway_[Name]_VEP.txt exists.
#
# Requirement: Rscript.exe is in the PATH, or adapt $RscriptPath below.
#
# Usage:
#   .\03b_run_phase2_all_pathways.ps1
# ==============================================================================

# FOCUS_BASE_DIR: project root containing the data folders (see README)
$BaseDir = if ($env:FOCUS_BASE_DIR) { $env:FOCUS_BASE_DIR } else { "path/to/FoCus" }
$RscriptPath = "C:/Program Files/R/R-4.5.2/bin/Rscript.exe"
$ScriptPath  = Join-Path $PSScriptRoot "03_pathway_pipeline.R"   # REPAIR
$PathwayListe = "$BaseDir/Pathway_Selection/Pathway_Liste.csv"
$LogDir = "$BaseDir/Repair/Batch_Logs"   # REPAIR

if (-not (Test-Path $LogDir)) {
    New-Item -ItemType Directory -Path $LogDir | Out-Null
}

$pathways = Import-Csv -Path $PathwayListe
if ($Only) { $pathways = $pathways | Where-Object { $Only -contains $_.pathway_name } }   # REPAIR: test-run filter

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

$summaryFile = Join-Path $LogDir "Phase2_Summary_REPAIR.csv"
$summary | Export-Csv -Path $summaryFile -NoTypeInformation -Encoding UTF8

Write-Host ""
Write-Host "Phase 2 completed for all $($pathways.Count) pathways." -ForegroundColor Yellow
Write-Host "Summary: $summaryFile"
