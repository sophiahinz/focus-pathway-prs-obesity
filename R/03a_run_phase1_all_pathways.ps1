# ==============================================================================
# Run_Phase1_AllPathways.ps1
#
# Runs the R script "03_pathway_pipeline.R"
# sequentially for ALL pathways from "Pathway_Liste.csv", each up to the
# breakpoint BEFORE the external VCF creation / VEP annotation (run_phase = 1).
#
# Requirement: Rscript.exe is in the PATH, or adapt $RscriptPath below.
#
# Usage:
#   .\03a_run_phase1_all_pathways.ps1
#
# After this run: perform the external VCF/VEP step for each pathway
# (see log note per pathway), then
# start 03b_run_phase2_all_pathways.ps1.
# ==============================================================================

# FOCUS_BASE_DIR: project root containing the data folders (see README)
$BaseDir = if ($env:FOCUS_BASE_DIR) { $env:FOCUS_BASE_DIR } else { "path/to/FoCus" }
$RscriptPath = "C:/Program Files/R/R-4.5.2/bin/Rscript.exe"
$ScriptPath  = Join-Path $PSScriptRoot "03_pathway_pipeline.R"
$PathwayListe = "$BaseDir/Pathway_Selection/Pathway_Liste.csv"
$LogDir = "$BaseDir/Batch_Logs"

if (-not (Test-Path $LogDir)) {
    New-Item -ItemType Directory -Path $LogDir | Out-Null
}

$pathways = Import-Csv -Path $PathwayListe

$summary = @()

foreach ($p in $pathways) {

    $name = $p.pathway_name
    $id   = $p.pathway_selection_id
    $logFile = Join-Path $LogDir "$($name)_Phase1.log"

    Write-Host "==> Phase 1: $name ($id)" -ForegroundColor Cyan

    & $RscriptPath $ScriptPath $name $id "1" *> $logFile
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

$summaryFile = Join-Path $LogDir "Phase1_Summary.csv"
$summary | Export-Csv -Path $summaryFile -NoTypeInformation -Encoding UTF8

Write-Host ""
Write-Host "Phase 1 completed for all $($pathways.Count) pathways." -ForegroundColor Yellow
Write-Host "Summary: $summaryFile"
Write-Host "Next step: external VCF creation + VEP annotation per pathway, then 03b_run_phase2_all_pathways.ps1"
