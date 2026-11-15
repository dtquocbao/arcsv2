#!/usr/bin/env bash
# Start ARC-SV Studio (U6). macOS / Ubuntu / Linux: ./scripts/start-arcsv.sh [--build]
# Top-level launchers: Start ARC-SV.command (macOS), Start ARC-SV.sh (Ubuntu/Linux).
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
chmod +x "$ROOT/scripts/start-arcsv.sh" "$ROOT/scripts/load-studio-image.sh" 2>/dev/null || true

BUILD=0
if [[ "${1:-}" == "--build" ]]; then
  BUILD=1
fi

STUDIO_IMAGE="${ARCSV2_STUDIO_IMAGE:-ghcr.io/arcsv2/arcsv2:studio}"
# Prefer ~/Documents (macOS + Ubuntu); fall back to ~/ARC-SV if Documents is missing.
if [[ -d "${HOME}/Documents" ]]; then
  DEFAULT_DATA="${HOME}/Documents/ARC-SV"
else
  DEFAULT_DATA="${HOME}/ARC-SV"
fi
PORT=8501
URL="http://127.0.0.1:${PORT}"

# Release image is linux/amd64. Apple Silicon and other hosts need emulation.
DOCKER_PLATFORM=()
HOST_ARCH="$(uname -m)"
case "$HOST_ARCH" in
  arm64|aarch64)
    DOCKER_PLATFORM=(--platform linux/amd64)
    ;;
esac
if [[ -n "${ARCSV2_DOCKER_PLATFORM:-}" ]]; then
  DOCKER_PLATFORM=(--platform "$ARCSV2_DOCKER_PLATFORM")
fi

if ! command -v docker >/dev/null 2>&1; then
  echo "Docker is not installed."
  echo "  macOS:  https://docs.docker.com/desktop/install/mac-install/"
  echo "  Ubuntu: https://docs.docker.com/engine/install/ubuntu/  (or Docker Desktop)"
  exit 1
fi
if ! docker info >/dev/null 2>&1; then
  echo "Docker is not running. Start Docker Desktop (macOS/Windows) or the docker service (Ubuntu)."
  exit 1
fi

mkdir -p "$DEFAULT_DATA"

pick_folder() {
  local default="$1"
  local chosen=""
  if [[ "$(uname -s)" == "Darwin" ]]; then
    chosen="$(osascript -e 'POSIX path of (choose folder with prompt "Choose where ARC-SV stores your data" default location POSIX file "'"$default"'")' 2>/dev/null || true)"
    if [[ -n "${chosen:-}" ]]; then
      echo "${chosen%/}"
      return
    fi
  fi
  if command -v zenity >/dev/null 2>&1; then
    chosen="$(zenity --file-selection --directory --title="Choose ARC-SV data folder" --filename="$default/" 2>/dev/null || true)"
    if [[ -n "${chosen:-}" ]]; then
      echo "$chosen"
      return
    fi
  fi
  if command -v kdialog >/dev/null 2>&1; then
    chosen="$(kdialog --getexistingdirectory "$default" 2>/dev/null || true)"
    if [[ -n "${chosen:-}" ]]; then
      echo "$chosen"
      return
    fi
  fi
  echo "$default"
}

DATA="${ARCSV2_DATA:-}"
if [[ -z "$DATA" ]]; then
  echo "Default data folder: $DEFAULT_DATA"
  DATA="$(pick_folder "$DEFAULT_DATA")"
fi
mkdir -p "$DATA"
export ARCSV2_DATA="$DATA"

if [[ "$BUILD" -eq 1 ]]; then
  echo "Building image from source (developer mode) …"
  "$ROOT/scripts/arcsv2.sh" build
  STUDIO_IMAGE="arcsv2:dev"
else
  HAVE_LOCAL=0
  for tag in arcsv2:studio "$STUDIO_IMAGE" arcsv2:dev; do
    if docker image inspect "$tag" >/dev/null 2>&1; then
      STUDIO_IMAGE="$tag"
      HAVE_LOCAL=1
      echo "Using local image: $STUDIO_IMAGE"
      break
    fi
  done
  if [[ "$HAVE_LOCAL" -eq 0 && -f "$ROOT/scripts/load-studio-image.sh" ]]; then
    echo "No local Studio image — trying load-studio-image.sh …"
    if bash "$ROOT/scripts/load-studio-image.sh"; then
      if docker image inspect arcsv2:studio >/dev/null 2>&1; then
        STUDIO_IMAGE="arcsv2:studio"
        HAVE_LOCAL=1
      fi
    fi
  fi
  if [[ "$HAVE_LOCAL" -eq 0 ]]; then
    echo "Pulling Studio image: $STUDIO_IMAGE"
    if ! docker pull "${DOCKER_PLATFORM[@]}" "$STUDIO_IMAGE"; then
      echo "Could not pull $STUDIO_IMAGE."
      echo "If you received a .tar image, run:"
      echo "  ./Load\\ Studio\\ Image.sh   (Ubuntu/Linux)"
      echo "  or double-click Load Studio Image.command (macOS)"
      echo "Developers: re-run with --build. Or set ARCSV2_STUDIO_IMAGE."
      exit 1
    fi
  fi
fi

if [[ ${#DOCKER_PLATFORM[@]} -gt 0 ]]; then
  echo "Note: host is $HOST_ARCH — running image as linux/amd64 (Docker Desktop uses emulation)."
fi

EXAMPLE_SRC="$ROOT/python/arcsv2/dashboard/demo/example_bam"
EXAMPLE_DST="$DATA/_bundled_example"
if [[ -f "$EXAMPLE_SRC/example.bam" && ! -f "$EXAMPLE_DST/example.bam" ]]; then
  echo "Copying bundled example BAM into data folder..."
  mkdir -p "$EXAMPLE_DST"
  cp -f "$EXAMPLE_SRC/"* "$EXAMPLE_DST/" 2>/dev/null || true
fi

echo "Starting ARC-SV Studio"
echo "  data:  $DATA"
echo "  image: $STUDIO_IMAGE"
echo "  url:   $URL"

(
  for _ in $(seq 1 90); do
    if curl -sf "$URL" >/dev/null 2>&1; then
      if command -v open >/dev/null 2>&1; then
        open "$URL"
      elif command -v xdg-open >/dev/null 2>&1; then
        xdg-open "$URL"
      fi
      exit 0
    fi
    sleep 2
  done
) &

EXTRA=()
if [[ "$BUILD" -eq 1 || -d "$ROOT/python/arcsv2" ]]; then
  EXTRA+=(-v "$ROOT:/work" -e PYTHONPATH=/work/python)
fi

docker run --rm -p "${PORT}:8501" \
  "${DOCKER_PLATFORM[@]}" \
  -v "$DATA:/data" \
  -e ARCSV2_RESULTS=/data \
  -e ARCSV2_CACHE=/data/ref_cache \
  -e ARCSV2_EXAMPLE_DIR=/data/_bundled_example \
  --name arcsv2-studio \
  "${EXTRA[@]}" \
  "$STUDIO_IMAGE" \
  dashboard --results /data --host 0.0.0.0
