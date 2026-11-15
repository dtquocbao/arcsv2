#!/bin/bash
# macOS: double-click to start ARC-SV Studio.
cd "$(dirname "$0")"
chmod +x ./scripts/start-arcsv.sh ./scripts/load-studio-image.sh \
  "./Load Studio Image.sh" "./Start ARC-SV.sh" 2>/dev/null || true
if [[ "${1:-}" == "--build" ]]; then
  exec ./scripts/start-arcsv.sh --build
fi
exec ./scripts/start-arcsv.sh
