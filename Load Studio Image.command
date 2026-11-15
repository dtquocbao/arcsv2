#!/bin/bash
# macOS: double-click to load the Docker image .tar before first start.
cd "$(dirname "$0")"
chmod +x ./scripts/load-studio-image.sh 2>/dev/null || true
./scripts/load-studio-image.sh
echo ""
echo "Press Enter to close…"
read -r _
