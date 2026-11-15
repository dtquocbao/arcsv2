# Build a Studio release package for sharing (Windows).
# Usage:
#   .\scripts\package_studio_release.ps1
#   .\scripts\package_studio_release.ps1 -IncludeImage   # also docker save (~GBs)
#   .\scripts\package_studio_release.ps1 -Version 0.1.0-studio
param(
    [string]$Version = "0.1.0-studio",
    [switch]$IncludeImage
)

$ErrorActionPreference = "Stop"
$Root = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")).Path
$Dist = Join-Path $Root "dist"
$OutDir = Join-Path $Dist "arcsv2-studio-$Version"
$Zip = Join-Path $Dist "arcsv2-studio-$Version.zip"
$ImageTar = Join-Path $Dist "arcsv2-studio-$Version-image.tar"
$StudioTag = "arcsv2:studio"
$GhcrTag = "ghcr.io/arcsv2/arcsv2:studio"

Write-Host "Packaging ARC-SV Studio $Version"
Write-Host "  root: $Root"

if (Test-Path -LiteralPath $OutDir) {
    Remove-Item -LiteralPath $OutDir -Recurse -Force
}
New-Item -ItemType Directory -Path (Join-Path $OutDir "scripts") -Force | Out-Null
New-Item -ItemType Directory -Path (Join-Path $OutDir "docs") -Force | Out-Null
New-Item -ItemType Directory -Path (Join-Path $OutDir "python\arcsv2\dashboard\demo\example_bam") -Force | Out-Null

$copies = @(
    @{ Src = "Start ARC-SV.bat"; Dst = "Start ARC-SV.bat" },
    @{ Src = "Start ARC-SV.command"; Dst = "Start ARC-SV.command" },
    @{ Src = "Start ARC-SV.sh"; Dst = "Start ARC-SV.sh" },
    @{ Src = "Load Studio Image.command"; Dst = "Load Studio Image.command" },
    @{ Src = "Load Studio Image.sh"; Dst = "Load Studio Image.sh" },
    @{ Src = "compose.studio.yaml"; Dst = "compose.studio.yaml" },
    @{ Src = "scripts\start-arcsv.ps1"; Dst = "scripts\start-arcsv.ps1" },
    @{ Src = "scripts\start-arcsv.sh"; Dst = "scripts\start-arcsv.sh" },
    @{ Src = "scripts\load-studio-image.ps1"; Dst = "scripts\load-studio-image.ps1" },
    @{ Src = "scripts\load-studio-image.sh"; Dst = "scripts\load-studio-image.sh" },
    @{ Src = "docs\QUESTIONS_FOR_DR_BO_SHORT.md"; Dst = "docs\QUESTIONS_FOR_DR_BO_SHORT.md" },
    @{ Src = "docs\TROUBLESHOOTING.md"; Dst = "docs\TROUBLESHOOTING.md" },
    @{ Src = "docs\CLUSTER.md"; Dst = "docs\CLUSTER.md" },
    @{ Src = "docs\OUTPUT.md"; Dst = "docs\OUTPUT.md" }
)

foreach ($item in $copies) {
    $src = Join-Path $Root $item.Src
    $dst = Join-Path $OutDir $item.Dst
    if (Test-Path -LiteralPath $src) {
        $parent = Split-Path -Parent $dst
        if (-not (Test-Path -LiteralPath $parent)) {
            New-Item -ItemType Directory -Path $parent -Force | Out-Null
        }
        Copy-Item -LiteralPath $src -Destination $dst -Force
    } else {
        Write-Host "  skip (missing): $($item.Src)" -ForegroundColor Yellow
    }
}

$exSrc = Join-Path $Root "python\arcsv2\dashboard\demo\example_bam"
if (Test-Path -LiteralPath $exSrc) {
    Copy-Item -Path (Join-Path $exSrc "*") -Destination (Join-Path $OutDir "python\arcsv2\dashboard\demo\example_bam") -Force
}

$forBo = @"
ARC-SV Studio — package for Dr. Bo
==================================
Version: $Version
Date: $(Get-Date -Format yyyy-MM-dd)

What this is
------------
A double-click browser app for structural variant calling (ARC-SV2).
No git clone or Python install is required.
Works on Windows, macOS, and Ubuntu (Docker required).

Install Docker
--------------
  Windows / macOS: Docker Desktop
    https://docs.docker.com/desktop/
  Ubuntu: Docker Engine (or Desktop)
    https://docs.docker.com/engine/install/ubuntu/
    Then: sudo usermod -aG docker `$USER   # log out/in once

What to do (5 minutes to first result)
--------------------------------------
1. Unzip this folder. Place the image .tar (if provided) next to it
   or inside it.

2. Load the Docker image once (skip if already loaded):
   Windows:  double-click  Load Studio Image.bat
   macOS:    double-click  Load Studio Image.command
   Ubuntu:   chmod +x "Load Studio Image.sh" "Start ARC-SV.sh"
             ./Load\ Studio\ Image.sh

3. Start Studio:
   Windows:  double-click  Start ARC-SV.bat
   macOS:    double-click  Start ARC-SV.command
             (if macOS blocks it: right-click → Open)
   Ubuntu:   ./Start\ ARC-SV.sh
             (or in Files: right-click → Allow Launching / Run as Program)

4. Accept the default data folder (Documents/ARC-SV) or pick another.
5. When the browser opens:
   - Setup  → confirm checks are green (optional: Download GRCh38 ~3 GB)
   - Try example → wait 2–5 minutes → open Results
   - Or New analysis → pick a BAM/CRAM → Quick test (chr20)

Apple Silicon Mac (M1/M2/M3/M4)
-------------------------------
The bundled image is linux/amd64. Docker Desktop runs it via emulation
automatically (the launcher passes --platform linux/amd64). First start
may be slower; that is expected.

Real data tip
-------------
Copy BAM/CRAM (+ index) into Documents/ARC-SV before New analysis.
For a first real run use Quick test (chr20). Whole genome needs more
disk/time and the GRCh38 reference (~3 GB, downloaded once).

Outputs
-------
Each run writes .tab and VCF (.vcf.gz + .tbi) under your data folder.
See docs/OUTPUT.md.

Cluster / 10k-sample cohort
---------------------------
See docs/CLUSTER.md (SLURM + Apptainer). The Studio app is for laptop
smoke tests; production batches use arcsv2 batch / Snakemake.

Open questions (paper vs original code)
---------------------------------------
docs/QUESTIONS_FOR_DR_BO_SHORT.md — short list of method questions
where the Cell 2024 Methods and the original ARC-SV code disagree.

Troubleshooting
---------------
docs/TROUBLESHOOTING.md and PLATFORM.txt
Common: Docker not running; pull failed → load the .tar image first;
corrupt GRCh38 cache → delete Documents/ARC-SV/ref_cache and re-download;
Ubuntu "Permission denied" on .sh → chmod +x the launcher scripts.

Contact
-------
Packaged from the ARC-SV2 usability (Studio) branch for review.
"@
Set-Content -LiteralPath (Join-Path $OutDir "FOR_DR_BO.txt") -Value $forBo -Encoding UTF8

$readme = @"
ARC-SV Studio $Version
====================
Windows · macOS · Ubuntu

1. Install Docker (Desktop on Windows/macOS; Engine or Desktop on Ubuntu).

2. If this package came with an image .tar, load it once:
   Windows:  Load Studio Image.bat
   macOS:    Load Studio Image.command
   Ubuntu:   chmod +x "Load Studio Image.sh" && ./Load\ Studio\ Image.sh

3. Start:
   Windows:  Start ARC-SV.bat
   macOS:    Start ARC-SV.command   (right-click → Open if Gatekeeper blocks)
   Ubuntu:   chmod +x "Start ARC-SV.sh" && ./Start\ ARC-SV.sh

4. Accept Documents/ARC-SV (or pick another folder).
5. Browser opens → Try example → Results.

See FOR_DR_BO.txt and PLATFORM.txt for details.
"@
Set-Content -LiteralPath (Join-Path $OutDir "README.txt") -Value $readme -Encoding UTF8

$platform = @"
Platform notes
==============

Windows
  - Docker Desktop required.
  - Start: Start ARC-SV.bat
  - Load image: Load Studio Image.bat

macOS (Intel or Apple Silicon)
  - Docker Desktop required.
  - Start: double-click Start ARC-SV.command
    If macOS says the app can't be opened: right-click → Open.
  - Load image: double-click Load Studio Image.command
  - Apple Silicon: image is linux/amd64; Docker emulates it (slightly slower).

Ubuntu / Linux
  - Install Docker Engine: https://docs.docker.com/engine/install/ubuntu/
    sudo apt-get update
    # follow Docker's install steps, then:
    sudo usermod -aG docker `$USER
    # log out and back in
  - Or install Docker Desktop for Linux.
  - After unzip (permissions may be lost in zip from Windows):
      chmod +x "Start ARC-SV.sh" "Load Studio Image.sh" scripts/*.sh
  - Start: ./Start\ ARC-SV.sh
  - Load image: ./Load\ Studio\ Image.sh
  - Optional GUI folder picker: sudo apt install zenity

Data folder
  Default: ~/Documents/ARC-SV  (or ~/ARC-SV if Documents is missing)
  Override: export ARCSV2_DATA=/path/to/folder

Image override
  export ARCSV2_STUDIO_IMAGE=arcsv2:studio
"@
Set-Content -LiteralPath (Join-Path $OutDir "PLATFORM.txt") -Value $platform -Encoding UTF8

$loadBat = @"
@echo off
cd /d "%~dp0"
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0scripts\load-studio-image.ps1"
if errorlevel 1 pause
"@
Set-Content -LiteralPath (Join-Path $OutDir "Load Studio Image.bat") -Value $loadBat -Encoding ASCII

# Prefer local studio tag; launchers also try loading sibling .tar
if ($IncludeImage) {
    Write-Host "Preparing Docker image…"
    $hasStudio = $false
    docker image inspect $StudioTag 2>$null | Out-Null
    if ($LASTEXITCODE -eq 0) { $hasStudio = $true }
    if (-not $hasStudio) {
        docker image inspect "arcsv2:dev" 2>$null | Out-Null
        if ($LASTEXITCODE -eq 0) {
            Write-Host "  Tagging arcsv2:dev → $StudioTag and $GhcrTag"
            docker tag arcsv2:dev $StudioTag
            docker tag arcsv2:dev $GhcrTag
            $hasStudio = $true
        }
    }
    if (-not $hasStudio) {
        throw "No arcsv2:studio or arcsv2:dev image found. Build with .\scripts\arcsv2.ps1 build first."
    }
    docker tag $StudioTag $GhcrTag 2>$null
    Write-Host "  Saving $StudioTag → $ImageTar (this can take several minutes)…"
    if (Test-Path -LiteralPath $ImageTar) { Remove-Item -LiteralPath $ImageTar -Force }
    docker save -o $ImageTar $StudioTag $GhcrTag
    $tarMb = [math]::Round((Get-Item -LiteralPath $ImageTar).Length / 1MB, 1)
    Write-Host "  Image tar: $tarMb MB"
    # Copy a small pointer file into the zip folder; keep the large tar beside the zip
    Set-Content -LiteralPath (Join-Path $OutDir "IMAGE.txt") -Value @"
Docker image for this release:
  $([IO.Path]::GetFileName($ImageTar))

Place that .tar next to this folder (or inside it), then load once:

  Windows:  Load Studio Image.bat
  macOS:    Load Studio Image.command
  Ubuntu:   ./Load\ Studio\ Image.sh

Then start Studio (Start ARC-SV.bat / .command / .sh).
The image is linux/amd64 and works on Windows, Intel Mac, Apple Silicon
(via Docker emulation), and Ubuntu x86_64 (and Ubuntu arm64 via emulation).
"@ -Encoding UTF8
}

# Point at an existing image tar even when -IncludeImage was not used this run.
if (-not (Test-Path -LiteralPath (Join-Path $OutDir "IMAGE.txt")) -and (Test-Path -LiteralPath $ImageTar)) {
    Set-Content -LiteralPath (Join-Path $OutDir "IMAGE.txt") -Value @"
Docker image for this release:
  $([IO.Path]::GetFileName($ImageTar))

Place that .tar next to this folder (or inside it), then load once:

  Windows:  Load Studio Image.bat
  macOS:    Load Studio Image.command
  Ubuntu:   ./Load\ Studio\ Image.sh

Then start Studio (Start ARC-SV.bat / .command / .sh).
"@ -Encoding UTF8
}

if (Test-Path -LiteralPath $Zip) { Remove-Item -LiteralPath $Zip -Force }
Write-Host "Creating zip $Zip …"
Compress-Archive -Path $OutDir -DestinationPath $Zip -Force

Write-Host ""
Write-Host "Done."
Write-Host "  Folder: $OutDir"
Write-Host "  Zip:    $Zip"
if (Test-Path -LiteralPath $ImageTar) {
    Write-Host "  Image:  $ImageTar"
    Write-Host ""
    Write-Host "Send to Dr. Bo:"
    Write-Host "  1) $Zip"
    Write-Host "  2) $ImageTar"
}
