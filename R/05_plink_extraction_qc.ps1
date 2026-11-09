# ==============================================================================
# 05_plink_extraction_qc.ps1
# Per-pathway SNP extraction and QC with PLINK 1.9, separately for the two
# FoCus recruitment arms (Ambu, MD), from the gene range files of
# 04_pathway_range_files.R (Supplementary Methods 1).
#
# Hinz et al. Pathway-restricted polygenic scores for energy-substrate
# metabolism and obesity (FoCus cohort). See README for inputs and run order.
# ==============================================================================
#
# For each file Working_Plink/pathway_<Name>_gene_ranges.txt and each cohort:
#   1) Extraction of all SNPs in the gene ranges of the pathway
#      (--extract range, --set-hh-missing)
#   2) QC, depending on the number of extracted SNPs:
#        > 5 SNPs : --maf 0.01 --geno 0.05 --mind 0.1
#        <= 5 SNPs: --mac 3    --geno 0.1  --mind 0.5  (small gene sets)
#
# Inputs (in <FOCUS_BASE_DIR>/Working_Plink):
#   FO_GSA_Ambu.bed/.bim/.fam, FO_GSA_MD.bed/.bim/.fam, pathway_*_gene_ranges.txt
# Outputs (in <FOCUS_BASE_DIR>/Working_Plink/plink_pathway_outputs):
#   Working_Plink_Pathway_<Name>_<Cohort>_raw.*   (extraction)
#   Working_Plink_Pathway_<Name>_<Cohort>_QC.*    (after QC; read by 07)
#   all_pathways_qc_summary.tsv                   (SNP/individual counts; read by 07)
#
# PLINK is taken from Working_Plink/plink.exe unless the environment variable
# PLINK points to the PLINK 1.9 executable.
# ==============================================================================

$BaseDir  = if ($env:FOCUS_BASE_DIR) { $env:FOCUS_BASE_DIR } else { "path/to/FoCus" }
$WorkDir  = Join-Path $BaseDir "Working_Plink"
$plink    = if ($env:PLINK) { $env:PLINK } else { Join-Path $WorkDir "plink.exe" }
$outdir   = Join-Path $WorkDir "plink_pathway_outputs"
$summary  = Join-Path $outdir "all_pathways_qc_summary.tsv"

New-Item -ItemType Directory -Force -Path $outdir | Out-Null

"PATHWAY`tCOHORT`tGENE_COUNT`tSNP_COUNT_RAW`tIND_COUNT_RAW`tQC_METHOD`tSNP_COUNT_QC`tIND_COUNT_QC" | Out-File $summary -Encoding ascii

$cohorts = @{
    "Ambu" = (Join-Path $WorkDir "FO_GSA_Ambu")
    "MD"   = (Join-Path $WorkDir "FO_GSA_MD")
}

$rangeFiles = Get-ChildItem -Path $WorkDir -Filter "pathway_*_gene_ranges.txt"
Write-Host "Found: $($rangeFiles.Count) pathway range files."

foreach ($rf in $rangeFiles) {

    # Pathway name from file name: pathway_XXX_gene_ranges.txt -> XXX
    $pathway  = $rf.BaseName -replace '^pathway_', '' -replace '_gene_ranges$', ''
    $geneList = $rf.FullName

    $geneCount = (
        Get-Content $geneList |
        Where-Object { $_ -notmatch "^\s*$" -and $_ -notmatch "^\s*CHR" }
    ).Count

    Write-Host "=== Pathway: $pathway ($geneCount genes) ==="

    foreach ($cohort in $cohorts.Keys) {

        $bfile = $cohorts[$cohort]

        $rawOut = Join-Path $outdir "Working_Plink_Pathway_$pathway`_$cohort`_raw"
        $qcOut  = Join-Path $outdir "Working_Plink_Pathway_$pathway`_$cohort`_QC"

        # 1) Extract all SNPs of the pathway's gene ranges
        & $plink `
            --bfile $bfile `
            --extract range $geneList `
            --set-hh-missing `
            --make-bed `
            --out $rawOut

        $snpCountRaw = if (Test-Path "$rawOut.bim") { (Get-Content "$rawOut.bim").Count } else { 0 }
        $indCountRaw = if (Test-Path "$rawOut.fam") { (Get-Content "$rawOut.fam").Count } else { 0 }

        # 2) QC: MAC-based thresholds for very small SNP sets, MAF-based otherwise
        if ($snpCountRaw -le 5) {

            $qcMethod = "MAC"

            & $plink `
                --bfile $rawOut `
                --set-hh-missing `
                --mac 3 `
                --geno 0.1 `
                --mind 0.5 `
                --make-bed `
                --out $qcOut

        } else {

            $qcMethod = "MAF"

            & $plink `
                --bfile $rawOut `
                --set-hh-missing `
                --maf 0.01 `
                --geno 0.05 `
                --mind 0.1 `
                --make-bed `
                --out $qcOut
        }

        $snpCountQc = if (Test-Path "$qcOut.bim") { (Get-Content "$qcOut.bim").Count } else { 0 }
        $indCountQc = if (Test-Path "$qcOut.fam") { (Get-Content "$qcOut.fam").Count } else { 0 }

        "$pathway`t$cohort`t$geneCount`t$snpCountRaw`t$indCountRaw`t$qcMethod`t$snpCountQc`t$indCountQc" | Add-Content $summary
    }
}

Write-Host "Done. Summary: $summary"
