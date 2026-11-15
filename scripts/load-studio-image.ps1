# Load ARC-SV Studio Docker image from a release .tar (offline-friendly).
param(
    [string]$TarPath = ""
)

$ErrorActionPreference = "Stop"
$Root = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")).Path
$StudioTag = "arcsv2:studio"
$GhcrTag = "ghcr.io/arcsv2/arcsv2:studio"

function Test-DockerRunning {
    if (-not (Get-Command docker -ErrorAction SilentlyContinue)) {
        Write-Host "Docker Desktop is not installed." -ForegroundColor Red
        Write-Host "Install from: https://docs.docker.com/desktop/" -ForegroundColor Yellow
        exit 1
    }
    cmd /c "docker info >nul 2>nul"
    if ($LASTEXITCODE -ne 0) {
        Write-Host "Docker Desktop is not running. Start it and try again." -ForegroundColor Red
        exit 1
    }
}

function Find-ImageTar {
    param([string]$Hint)
    if ($Hint -and (Test-Path -LiteralPath $Hint)) {
        return (Resolve-Path -LiteralPath $Hint).Path
    }
    $candidates = @()
    $candidates += Get-ChildItem -LiteralPath $Root -Filter "arcsv2-studio*-image.tar" -ErrorAction SilentlyContinue
    $parent = Split-Path -Parent $Root
    if ($parent) {
        $candidates += Get-ChildItem -LiteralPath $parent -Filter "arcsv2-studio*-image.tar" -ErrorAction SilentlyContinue
    }
    $candidates += Get-ChildItem -LiteralPath $Root -Filter "*.tar" -ErrorAction SilentlyContinue |
        Where-Object { $_.Name -match "arcsv2|studio" }
    $hit = $candidates | Sort-Object LastWriteTime -Descending | Select-Object -First 1
    if ($hit) { return $hit.FullName }
    return $null
}

Test-DockerRunning

docker image inspect $StudioTag 2>$null | Out-Null
if ($LASTEXITCODE -eq 0) {
    Write-Host "Image already present: $StudioTag"
    docker tag $StudioTag $GhcrTag 2>$null
    Write-Host "You can Start ARC-SV.bat now."
    exit 0
}

$tar = Find-ImageTar -Hint $TarPath
if (-not $tar) {
    Write-Host "No image .tar found next to this package." -ForegroundColor Red
    Write-Host "Expected something like: arcsv2-studio-*-image.tar" -ForegroundColor Yellow
    Write-Host "Place it in:" -ForegroundColor Yellow
    Write-Host "  $Root"
    Write-Host "or the parent folder, then re-run this script."
    exit 1
}

Write-Host "Loading Docker image from:"
Write-Host "  $tar"
Write-Host "(this can take several minutes)…"
docker load -i $tar
if ($LASTEXITCODE -ne 0) {
    Write-Host "docker load failed." -ForegroundColor Red
    exit 1
}

docker image inspect $StudioTag 2>$null | Out-Null
if ($LASTEXITCODE -ne 0) {
    docker image inspect $GhcrTag 2>$null | Out-Null
    if ($LASTEXITCODE -eq 0) {
        docker tag $GhcrTag $StudioTag
    } else {
        docker image inspect "arcsv2:dev" 2>$null | Out-Null
        if ($LASTEXITCODE -eq 0) {
            docker tag arcsv2:dev $StudioTag
            docker tag arcsv2:dev $GhcrTag
        }
    }
}
docker tag $StudioTag $GhcrTag 2>$null

Write-Host ""
Write-Host "Loaded. Next: double-click Start ARC-SV.bat"
