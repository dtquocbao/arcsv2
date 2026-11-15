#!/usr/bin/env bash
# Ubuntu / Linux launcher (double-click in Files if "Run as program" is enabled,
# or: chmod +x "Start ARC-SV.sh" && ./Start\ ARC-SV.sh)
cd "$(dirname "$0")"
chmod +x ./scripts/start-arcsv.sh ./scripts/load-studio-image.sh 2>/dev/null || true
if [[ "${1:-}" == "--build" ]]; then
  exec ./scripts/start-arcsv.sh --build
fi
exec ./scripts/start-arcsv.sh
