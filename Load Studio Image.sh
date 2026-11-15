#!/usr/bin/env bash
# Ubuntu / Linux: load the Docker image .tar before first start.
cd "$(dirname "$0")"
chmod +x ./scripts/load-studio-image.sh 2>/dev/null || true
exec ./scripts/load-studio-image.sh "$@"
