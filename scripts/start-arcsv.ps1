# Start ARC-SV Studio (U6). Double-click "Start ARC-SV.bat".
# Default: pull prebuilt ghcr.io/arcsv2/arcsv2:studio
# Developer:  .\scripts\start-arcsv.ps1 -Build
param(
    [switch]$Build
)

$ErrorActionPreference = "Stop"
$Root = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")).Path
Set-Location -LiteralPath $Root

$StudioImage = if ($env:ARCSV2_STUDIO_IMAGE) { $env:ARCSV2_STUDIO_IMAGE } else { "ghcr.io/arcsv2/arcsv2:studio" }
$DefaultData = Join-Path ([Environment]::GetFolderPath("MyDocuments")) "ARC-SV"
$Port = 8501
$Url = "http://127.0.0.1:$Port"

function Test-DockerRunning {
    if (-not (Get-Command docker -ErrorAction SilentlyContinue)) {
        Write-Host ""
        Write-Host "Docker Desktop is not installed." -ForegroundColor Red
        Write-Host "Install from: https://docs.docker.com/desktop/" -ForegroundColor Yellow
        Read-Host "Press Enter to close"
        exit 1
    }
    # Docker prints engine warnings to stderr; do not treat those as terminating errors.
    cmd /c "docker info >nul 2>nul"
    if ($LASTEXITCODE -ne 0) {
        Write-Host ""
        Write-Host "Docker Desktop is installed but not running." -ForegroundColor Red
        Write-Host "Start Docker Desktop, wait until it is ready, then run this again." -ForegroundColor Yellow
        Write-Host "Install: https://docs.docker.com/desktop/" -ForegroundColor Yellow
        Read-Host "Press Enter to close"
        exit 1
    }
}

function Choose-DataFolder {
    param([string]$Default)
    Add-Type -AssemblyName System.Windows.Forms | Out-Null
    $dialog = New-Object System.Windows.Forms.FolderBrowserDialog
    $dialog.Description = "Choose where ARC-SV stores your data (BAM/CRAM and results)"
    $dialog.SelectedPath = $Default
    $dialog.ShowNewFolderButton = $true
    $result = $dialog.ShowDialog()
    if ($result -eq [System.Windows.Forms.DialogResult]::OK -and $dialog.SelectedPath) {
        return $dialog.SelectedPath
    }
    return $Default
}

Test-DockerRunning

if (-not (Test-Path -LiteralPath $DefaultData)) {
    New-Item -ItemType Directory -Path $DefaultData -Force | Out-Null
}

$DataFolder = $env:ARCSV2_DATA
if ([string]::IsNullOrWhiteSpace($DataFolder)) {
    Write-Host "Default data folder: $DefaultData"
    Write-Host "A folder picker will open - choose this folder or another."
    $DataFolder = Choose-DataFolder -Default $DefaultData
}
if (-not (Test-Path -LiteralPath $DataFolder)) {
    New-Item -ItemType Directory -Path $DataFolder -Force | Out-Null
}
$DataFull = (Resolve-Path -LiteralPath $DataFolder).Path
$env:ARCSV2_DATA = $DataFull

if ($Build) {
    Write-Host "Building image from source (developer mode)..."
    & "$Root\scripts\arcsv2.ps1" build
    $StudioImage = "arcsv2:dev"
} else {
    $haveLocal = $false
    foreach ($tag in @("arcsv2:studio", $StudioImage, "arcsv2:dev")) {
        docker image inspect $tag 2>$null | Out-Null
        if ($LASTEXITCODE -eq 0) {
            $StudioImage = $tag
            $haveLocal = $true
            Write-Host "Using local image: $StudioImage"
            break
        }
    }
    if (-not $haveLocal) {
        $loadScript = Join-Path $PSScriptRoot "load-studio-image.ps1"
        if (Test-Path -LiteralPath $loadScript) {
            Write-Host "No local Studio image — trying load-studio-image.ps1 …"
            & $loadScript
            if ($LASTEXITCODE -eq 0) {
                docker image inspect "arcsv2:studio" 2>$null | Out-Null
                if ($LASTEXITCODE -eq 0) {
                    $StudioImage = "arcsv2:studio"
                    $haveLocal = $true
                }
            }
        }
    }
    if (-not $haveLocal) {
        Write-Host "Pulling Studio image: $StudioImage"
        cmd /c "docker pull $StudioImage"
        if ($LASTEXITCODE -ne 0) {
            Write-Host ""
            Write-Host "Could not pull $StudioImage." -ForegroundColor Red
            Write-Host "If you received a .tar image, run Load Studio Image.bat first." -ForegroundColor Yellow
            Write-Host "Developers: re-run with -Build to build from this checkout." -ForegroundColor Yellow
            Write-Host "Or set ARCSV2_STUDIO_IMAGE to a reachable tag." -ForegroundColor Yellow
            Read-Host "Press Enter to close"
            exit 1
        }
    }
}

# Seed the tiny example BAM into the data folder so Studio works even when the
# image was built before package-data included demo/example_bam.
$ExampleSrc = Join-Path $Root "python\arcsv2\dashboard\demo\example_bam"
$ExampleDst = Join-Path $DataFull "_bundled_example"
if ((Test-Path -LiteralPath (Join-Path $ExampleSrc "example.bam")) -and
    -not (Test-Path -LiteralPath (Join-Path $ExampleDst "example.bam"))) {
    Write-Host "Copying bundled example BAM into data folder..."
    New-Item -ItemType Directory -Path $ExampleDst -Force | Out-Null
    Copy-Item -Path (Join-Path $ExampleSrc "*") -Destination $ExampleDst -Force
}

Write-Host "Starting ARC-SV Studio"
Write-Host "  data:  $DataFull"
Write-Host "  image: $StudioImage"
Write-Host "  url:   $Url"

$job = Start-Job -ScriptBlock {
    param($u)
    $deadline = (Get-Date).AddSeconds(180)
    while ((Get-Date) -lt $deadline) {
        try {
            $r = Invoke-WebRequest -Uri $u -UseBasicParsing -TimeoutSec 2
            if ($r.StatusCode -ge 200) {
                Start-Process $u
                return
            }
        } catch {
            Start-Sleep -Seconds 2
        }
    }
} -ArgumentList $Url

# Argument array avoids backtick line-continuation bugs (trailing spaces / smart quotes).
$Volume = $DataFull + ':/data'
$Publish = [string]$Port + ':8501'
$DockerArgs = @(
    'run', '--rm',
    '-p', $Publish,
    '-v', $Volume,
    '-e', 'ARCSV2_RESULTS=/data',
    '-e', 'ARCSV2_CACHE=/data/ref_cache',
    '-e', 'ARCSV2_EXAMPLE_DIR=/data/_bundled_example',
    '--name', 'arcsv2-studio'
)
# Developer checkout: prefer live python tree (includes example BAM + Studio fixes).
if ($Build -or (Test-Path -LiteralPath (Join-Path $Root "python\arcsv2"))) {
    $DockerArgs += @('-v', ($Root + ':/work'), '-e', 'PYTHONPATH=/work/python')
}
$DockerArgs += @(
    $StudioImage,
    'dashboard', '--results', '/data', '--host', '0.0.0.0'
)
& docker @DockerArgs
$exitCode = $LASTEXITCODE

Receive-Job $job -ErrorAction SilentlyContinue | Out-Null
Remove-Job $job -Force -ErrorAction SilentlyContinue
exit $exitCode
