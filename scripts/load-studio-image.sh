#!/usr/bin/env bash
# Load ARC-SV Studio Docker image from a release .tar (offline-friendly).
# macOS: Load Studio Image.command · Ubuntu/Linux: Load Studio Image.sh
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STUDIO_TAG="arcsv2:studio"
GHCR_TAG="ghcr.io/arcsv2/arcsv2:studio"
HINT="${1:-}"

if ! command -v docker >/dev/null 2>&1; then
  echo "Docker is not installed."
  echo "  macOS:  https://docs.docker.com/desktop/install/mac-install/"
  echo "  Ubuntu: https://docs.docker.com/engine/install/ubuntu/"
  exit 1
fi
if ! docker info >/dev/null 2>&1; then
  echo "Docker is not running. Start Docker Desktop (macOS) or: sudo systemctl start docker (Ubuntu)."
  exit 1
fi

if docker image inspect "$STUDIO_TAG" >/dev/null 2>&1; then
  echo "Image already present: $STUDIO_TAG"
  docker tag "$STUDIO_TAG" "$GHCR_TAG" 2>/dev/null || true
  echo "You can start Studio now."
  exit 0
fi

find_tar() {
  if [[ -n "$HINT" && -f "$HINT" ]]; then
    echo "$HINT"
    return
  fi
  local f
  f="$(ls -t "$ROOT"/arcsv2-studio*-image.tar 2>/dev/null | head -1 || true)"
  if [[ -n "$f" ]]; then echo "$f"; return; fi
  f="$(ls -t "$(dirname "$ROOT")"/arcsv2-studio*-image.tar 2>/dev/null | head -1 || true)"
  if [[ -n "$f" ]]; then echo "$f"; return; fi
  f="$(ls -t "$ROOT"/*studio*.tar "$ROOT"/arcsv2*.tar 2>/dev/null | head -1 || true)"
  if [[ -n "$f" ]]; then echo "$f"; return; fi
  echo ""
}

TAR="$(find_tar)"
if [[ -z "$TAR" ]]; then
  echo "No image .tar found next to this package."
  echo "Expected: arcsv2-studio-*-image.tar in $ROOT or its parent."
  exit 1
fi

echo "Loading Docker image from: $TAR"
docker load -i "$TAR"

if ! docker image inspect "$STUDIO_TAG" >/dev/null 2>&1; then
  if docker image inspect "$GHCR_TAG" >/dev/null 2>&1; then
    docker tag "$GHCR_TAG" "$STUDIO_TAG"
  elif docker image inspect arcsv2:dev >/dev/null 2>&1; then
    docker tag arcsv2:dev "$STUDIO_TAG"
    docker tag arcsv2:dev "$GHCR_TAG"
  fi
fi
docker tag "$STUDIO_TAG" "$GHCR_TAG" 2>/dev/null || true
echo "Loaded. Next: ./Start ARC-SV.command or ./scripts/start-arcsv.sh"
