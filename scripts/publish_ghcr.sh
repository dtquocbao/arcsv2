#!/usr/bin/env bash
# Publish frozen arcsv2:ml0 (caller-v1.6.1) to GitHub Container Registry.
#
# Usage (laptop with Docker Desktop running):
#   export GHCR_USER=dtquocbao          # GitHub username / org that owns the package
#   export GHCR_TOKEN=ghp_...           # PAT with write:packages (and read:packages)
#   # or: echo "$GHCR_TOKEN" | docker login ghcr.io -u "$GHCR_USER" --password-stdin
#   bash scripts/publish_ghcr.sh
#
# Tags pushed:
#   ghcr.io/$GHCR_USER/arcsv2:1.6.1
#   ghcr.io/$GHCR_USER/arcsv2:latest-dev
#   ghcr.io/$GHCR_USER/arcsv2:ml0
#
# After push, make the package public (one-time) so Grace can pull without a token:
#   gh api --method PUT -H "Accept: application/vnd.github+json" \
#     "/user/packages/container/arcsv2/visibility" -f visibility=public
# Or: GitHub → Packages → arcsv2 → Package settings → Change visibility → Public
set -euo pipefail

GHCR_USER="${GHCR_USER:-dtquocbao}"
SRC_IMAGE="${SRC_IMAGE:-arcsv2:ml0}"
VERSION_TAG="${VERSION_TAG:-1.6.1}"
REGISTRY="ghcr.io/${GHCR_USER}/arcsv2"

if ! command -v docker >/dev/null 2>&1; then
  echo "docker not found" >&2
  exit 1
fi
if ! docker info >/dev/null 2>&1; then
  echo "Docker daemon is not running. Start Docker Desktop and retry." >&2
  exit 1
fi
if ! docker image inspect "$SRC_IMAGE" >/dev/null 2>&1; then
  echo "Source image missing: $SRC_IMAGE" >&2
  echo "Retag from current freeze, e.g.: docker tag arcsv2:dev arcsv2:ml0" >&2
  exit 1
fi

if [[ -n "${GHCR_TOKEN:-}" ]]; then
  echo "$GHCR_TOKEN" | docker login ghcr.io -u "$GHCR_USER" --password-stdin
elif ! docker pull "ghcr.io/${GHCR_USER}/hello-world:nonexistent" >/dev/null 2>&1; then
  # Probe login: if not logged in, docker info about auth may still pass; try a no-op login check.
  if ! grep -q "ghcr.io" "${HOME}/.docker/config.json" 2>/dev/null \
     && ! grep -q "ghcr.io" "${USERPROFILE:-}/.docker/config.json" 2>/dev/null; then
    echo "Not logged in to ghcr.io. Set GHCR_TOKEN or run: docker login ghcr.io -u $GHCR_USER" >&2
    exit 1
  fi
fi

echo "Tagging $SRC_IMAGE → ${REGISTRY}:{${VERSION_TAG},latest-dev,ml0}"
docker tag "$SRC_IMAGE" "${REGISTRY}:${VERSION_TAG}"
docker tag "$SRC_IMAGE" "${REGISTRY}:latest-dev"
docker tag "$SRC_IMAGE" "${REGISTRY}:ml0"

echo "Pushing…"
docker push "${REGISTRY}:${VERSION_TAG}"
docker push "${REGISTRY}:latest-dev"
docker push "${REGISTRY}:ml0"

echo ""
echo "Pushed:"
echo "  ${REGISTRY}:${VERSION_TAG}"
echo "  ${REGISTRY}:latest-dev"
echo "  ${REGISTRY}:ml0"
echo ""
echo "Make package public (so Grace singularity pull needs no token):"
echo "  gh api --method PUT -H 'Accept: application/vnd.github+json' \\"
echo "    /user/packages/container/arcsv2/visibility -f visibility=public"
echo "Or open: https://github.com/users/${GHCR_USER}/packages/container/package/arcsv2"
