# Publish frozen arcsv2:ml0 (caller-v1.6.1) to GHCR.
# Usage:
#   $env:GHCR_USER = "dtquocbao"
#   $env:GHCR_TOKEN = "ghp_..."   # PAT: write:packages, read:packages
#   .\scripts\publish_ghcr.ps1
param(
    [string]$GhcrUser = $(if ($env:GHCR_USER) { $env:GHCR_USER } else { "dtquocbao" }),
    [string]$SrcImage = $(if ($env:SRC_IMAGE) { $env:SRC_IMAGE } else { "arcsv2:ml0" }),
    [string]$VersionTag = $(if ($env:VERSION_TAG) { $env:VERSION_TAG } else { "1.6.1" })
)

$ErrorActionPreference = "Stop"
$Registry = "ghcr.io/$GhcrUser/arcsv2"

docker info 2>$null | Out-Null
if ($LASTEXITCODE -ne 0) {
    throw "Docker daemon is not running. Start Docker Desktop and retry."
}
docker image inspect $SrcImage 2>$null | Out-Null
if ($LASTEXITCODE -ne 0) {
    throw "Source image missing: $SrcImage. Retag with: docker tag arcsv2:dev arcsv2:ml0"
}

if ($env:GHCR_TOKEN) {
    $env:GHCR_TOKEN | docker login ghcr.io -u $GhcrUser --password-stdin
    if ($LASTEXITCODE -ne 0) { throw "docker login ghcr.io failed" }
}

Write-Host "Tagging $SrcImage → ${Registry}:{${VersionTag},latest-dev,ml0}"
docker tag $SrcImage "${Registry}:${VersionTag}"
docker tag $SrcImage "${Registry}:latest-dev"
docker tag $SrcImage "${Registry}:ml0"

Write-Host "Pushing…"
docker push "${Registry}:${VersionTag}"
if ($LASTEXITCODE -ne 0) { throw "push :$VersionTag failed" }
docker push "${Registry}:latest-dev"
if ($LASTEXITCODE -ne 0) { throw "push :latest-dev failed" }
docker push "${Registry}:ml0"
if ($LASTEXITCODE -ne 0) { throw "push :ml0 failed" }

Write-Host ""
Write-Host "Pushed:"
Write-Host "  ${Registry}:${VersionTag}"
Write-Host "  ${Registry}:latest-dev"
Write-Host "  ${Registry}:ml0"
Write-Host ""
Write-Host "Make package public (Grace pull without token):"
Write-Host "  gh api --method PUT -H `"Accept: application/vnd.github+json`" ``"
Write-Host "    /user/packages/container/arcsv2/visibility -f visibility=public"
Write-Host "Or: https://github.com/users/$GhcrUser/packages/container/package/arcsv2"
