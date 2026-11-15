#!/usr/bin/env bash
# Upload dist/release-staging assets to GitHub Release benchmark-2026-11.
#   gh auth login && bash scripts/upload_benchmark_release.sh
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPO="${REPO:-dtquocbao/arcsv2}"
TAG="${TAG:-benchmark-2026-11}"
STAGING="${STAGING:-$ROOT/dist/release-staging}"

gh auth status >/dev/null

files=(
  "$STAGING/SHA256SUMS"
  "$STAGING/dataset.tsv"
  "$STAGING/dataset_v0.parquet"
  "$STAGING/dataset_v0.tsv.gz"
  "$STAGING/ml0_merged.tar.gz"
  "$STAGING/export_history.tar.gz"
  "$STAGING/s4_host.log.gz"
  "$STAGING/s4_snakemake.log.gz"
)
for f in "${files[@]}"; do
  [[ -f "$f" ]] || { echo "Missing $f" >&2; exit 1; }
done

NOTES="Benchmark evidence and large ML-0 assets for ARC-SV2 (caller-v1.6.1).

See docs/DATA.md and docs/ORIGINAL_COMPARISON.md.

Verify:
\`\`\`
gh release download ${TAG} -R ${REPO} -D /tmp/arcsv2-data
cd /tmp/arcsv2-data && sha256sum -c SHA256SUMS
\`\`\`"

if gh release view "$TAG" --repo "$REPO" >/dev/null 2>&1; then
  echo "Release $TAG exists — uploading/replacing assets…"
  gh release upload "$TAG" --repo "$REPO" --clobber "${files[@]}"
else
  echo "Creating release $TAG…"
  gh release create "$TAG" --repo "$REPO" --title "$TAG" --notes "$NOTES" "${files[@]}"
fi

echo "Done: https://github.com/${REPO}/releases/tag/${TAG}"
