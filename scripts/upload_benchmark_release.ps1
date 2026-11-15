# Upload dist/release-staging assets to GitHub Release benchmark-2026-11.
# Requires: gh auth login (repo scope). You can push separately.
#
#   .\scripts\upload_benchmark_release.ps1
param(
    [string]$Repo = "dtquocbao/arcsv2",
    [string]$Tag = "benchmark-2026-11",
    [string]$Staging = ""
)

$ErrorActionPreference = "Stop"
$Root = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")).Path
if (-not $Staging) {
    $Staging = Join-Path $Root "dist\release-staging"
}

gh auth status 2>$null | Out-Null
if ($LASTEXITCODE -ne 0) {
    Write-Host "Run: gh auth login" -ForegroundColor Yellow
    exit 1
}

$required = @(
    "SHA256SUMS",
    "dataset.tsv",
    "dataset_v0.parquet",
    "dataset_v0.tsv.gz",
    "ml0_merged.tar.gz",
    "export_history.tar.gz",
    "s4_host.log.gz",
    "s4_snakemake.log.gz"
)
$files = @()
foreach ($name in $required) {
    $path = Join-Path $Staging $name
    if (-not (Test-Path -LiteralPath $path)) {
        throw "Missing staging file: $path"
    }
    $files += $path
}

$notes = @"
Benchmark evidence and large ML-0 assets for ARC-SV2 (``caller-v1.6.1``).

See ``docs/DATA.md`` and ``docs/ORIGINAL_COMPARISON.md``.

Verify:
``````
gh release download $Tag -R $Repo -D /tmp/arcsv2-data
cd /tmp/arcsv2-data && sha256sum -c SHA256SUMS
``````
"@

$exists = $false
gh release view $Tag --repo $Repo 2>$null | Out-Null
if ($LASTEXITCODE -eq 0) { $exists = $true }

if ($exists) {
    Write-Host "Release $Tag exists — uploading/replacing assets…"
    gh release upload $Tag --repo $Repo --clobber @files
} else {
    Write-Host "Creating release $Tag…"
    gh release create $Tag --repo $Repo --title $Tag --notes $notes @files
}

Write-Host ""
Write-Host "Done: https://github.com/$Repo/releases/tag/$Tag"
