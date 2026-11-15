#!/usr/bin/env bash
# Build a Studio release zip (and optional docker image tar) for sharing.
# Works when packaging on macOS/Linux; Windows maintainers prefer the .ps1.
# Usage:
#   bash scripts/package_studio_release.sh [VERSION]
#   INCLUDE_IMAGE=1 bash scripts/package_studio_release.sh 0.1.0-studio
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VER="${1:-0.1.0-studio}"
OUT_DIR="${ROOT}/dist/arcsv2-studio-${VER}"
ZIP="${ROOT}/dist/arcsv2-studio-${VER}.zip"
IMAGE_TAR="${ROOT}/dist/arcsv2-studio-${VER}-image.tar"
INCLUDE_IMAGE="${INCLUDE_IMAGE:-0}"
STUDIO_TAG="arcsv2:studio"
GHCR_TAG="ghcr.io/arcsv2/arcsv2:studio"
DATE="$(date +%Y-%m-%d)"

rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR/scripts" "$OUT_DIR/docs" \
  "$OUT_DIR/python/arcsv2/dashboard/demo/example_bam"

cp "$ROOT/Start ARC-SV.bat" "$OUT_DIR/"
cp "$ROOT/Start ARC-SV.command" "$OUT_DIR/"
cp "$ROOT/Start ARC-SV.sh" "$OUT_DIR/"
cp "$ROOT/Load Studio Image.command" "$OUT_DIR/"
cp "$ROOT/Load Studio Image.sh" "$OUT_DIR/"
cp "$ROOT/compose.studio.yaml" "$OUT_DIR/"
cp "$ROOT/scripts/start-arcsv.ps1" "$OUT_DIR/scripts/"
cp "$ROOT/scripts/start-arcsv.sh" "$OUT_DIR/scripts/"
cp "$ROOT/scripts/load-studio-image.ps1" "$OUT_DIR/scripts/"
cp "$ROOT/scripts/load-studio-image.sh" "$OUT_DIR/scripts/"
chmod +x \
  "$OUT_DIR/scripts/start-arcsv.sh" \
  "$OUT_DIR/scripts/load-studio-image.sh" \
  "$OUT_DIR/Start ARC-SV.command" \
  "$OUT_DIR/Start ARC-SV.sh" \
  "$OUT_DIR/Load Studio Image.command" \
  "$OUT_DIR/Load Studio Image.sh" || true

for doc in QUESTIONS_FOR_DR_BO_SHORT.md TROUBLESHOOTING.md CLUSTER.md OUTPUT.md; do
  if [[ -f "$ROOT/docs/$doc" ]]; then
    cp "$ROOT/docs/$doc" "$OUT_DIR/docs/"
  fi
done

if [[ -d "$ROOT/python/arcsv2/dashboard/demo/example_bam" ]]; then
  cp -f "$ROOT/python/arcsv2/dashboard/demo/example_bam/"* \
    "$OUT_DIR/python/arcsv2/dashboard/demo/example_bam/" || true
fi

cat > "$OUT_DIR/Load Studio Image.bat" <<'EOF'
@echo off
cd /d "%~dp0"
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0scripts\load-studio-image.ps1"
if errorlevel 1 pause
EOF

cat > "$OUT_DIR/FOR_DR_BO.txt" <<EOF
ARC-SV Studio — package for Dr. Bo
==================================
Version: ${VER}
Date: ${DATE}

Works on Windows, macOS, and Ubuntu (Docker required).

Install Docker
--------------
  Windows / macOS: https://docs.docker.com/desktop/
  Ubuntu: https://docs.docker.com/engine/install/ubuntu/
    Then: sudo usermod -aG docker \$USER   # log out/in once

Load image once (if you received a .tar)
----------------------------------------
  Windows:  Load Studio Image.bat
  macOS:    Load Studio Image.command
  Ubuntu:   chmod +x "Load Studio Image.sh" "Start ARC-SV.sh"
            ./Load\\ Studio\\ Image.sh

Start Studio
------------
  Windows:  Start ARC-SV.bat
  macOS:    Start ARC-SV.command  (right-click → Open if blocked)
  Ubuntu:   ./Start\\ ARC-SV.sh

Then: Try example → Results. See PLATFORM.txt for Apple Silicon / Ubuntu notes.
Questions: docs/QUESTIONS_FOR_DR_BO_SHORT.md
EOF

cat > "$OUT_DIR/README.txt" <<EOF
ARC-SV Studio ${VER}
====================
Windows · macOS · Ubuntu

1. Install Docker.
2. Load image .tar once (Load Studio Image.bat / .command / .sh).
3. Start (Start ARC-SV.bat / .command / .sh).
4. Browser → Try example → Results.

See FOR_DR_BO.txt and PLATFORM.txt.
EOF

cat > "$OUT_DIR/PLATFORM.txt" <<'EOF'
Platform notes
==============

Windows:  Start ARC-SV.bat · Load Studio Image.bat · Docker Desktop
macOS:    Start ARC-SV.command · Load Studio Image.command · Docker Desktop
          Apple Silicon uses linux/amd64 emulation automatically.
Ubuntu:   chmod +x "Start ARC-SV.sh" "Load Studio Image.sh" scripts/*.sh
          ./Start\ ARC-SV.sh · ./Load\ Studio\ Image.sh
          Docker Engine or Desktop; add user to docker group.
EOF

if [[ "$INCLUDE_IMAGE" == "1" ]]; then
  if ! docker image inspect "$STUDIO_TAG" >/dev/null 2>&1; then
    if docker image inspect arcsv2:dev >/dev/null 2>&1; then
      docker tag arcsv2:dev "$STUDIO_TAG"
      docker tag arcsv2:dev "$GHCR_TAG"
    else
      echo "No arcsv2:studio or arcsv2:dev image. Build first." >&2
      exit 1
    fi
  fi
  docker tag "$STUDIO_TAG" "$GHCR_TAG" 2>/dev/null || true
  echo "Saving Docker image to $IMAGE_TAR …"
  rm -f "$IMAGE_TAR"
  docker save -o "$IMAGE_TAR" "$STUDIO_TAG" "$GHCR_TAG"
  cat > "$OUT_DIR/IMAGE.txt" <<EOF
Docker image: $(basename "$IMAGE_TAR")
Load with Load Studio Image.bat / .command / .sh before starting Studio.
linux/amd64 — works on Windows, Intel/Apple Silicon Mac, Ubuntu (emulation on ARM).
EOF
fi

(cd "$(dirname "$OUT_DIR")" && rm -f "$(basename "$ZIP")" && zip -r "$(basename "$ZIP")" "$(basename "$OUT_DIR")")
echo "Wrote $ZIP"
if [[ "$INCLUDE_IMAGE" == "1" ]]; then
  echo "Wrote $IMAGE_TAR"
fi
