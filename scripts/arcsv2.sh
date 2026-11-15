#!/usr/bin/env bash
# Run arcsv2 via Compose (project name arcsv2). Repo is /work, data is /data.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
compose=(docker compose -f "$root/compose.yaml" --project-directory "$root")

ensure_data() {
  # Named volumes are created root-owned; mambauser must be able to write /data.
  "${compose[@]}" --profile init run --rm init-data >/dev/null
}

# Map a host path under the repo to the container mount /work/...
to_container_path() {
  local path="$1"
  case "$path" in
    /*) printf '%s' "$path"; return ;;
  esac
  local full
  full="$(cd "$(dirname "$path")" 2>/dev/null && pwd)/$(basename "$path")" || full="$path"
  case "$full" in
    "$root"|"$root"/*)
      local rel="${full#"$root"}"
      rel="${rel#/}"
      if [[ -z "$rel" ]]; then
        printf '/work'
      else
        printf '/work/%s' "$rel"
      fi
      ;;
    *) printf '%s' "$path" ;;
  esac
}

dashboard_args() {
  # stdout: words for the container command (starts with dashboard)
  local -a out=(dashboard)
  local has_host=0 has_results=0
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --host) has_host=1; out+=("$1"); shift; [[ $# -gt 0 ]] && out+=("$1") && shift ;;
      --results)
        has_results=1
        out+=("$1")
        shift
        if [[ $# -gt 0 ]]; then
          out+=("$(to_container_path "$1")")
          shift
        fi
        ;;
      *) out+=("$1"); shift ;;
    esac
  done
  # Fetch/Snakemake write under the named volume, not the repo bind.
  [[ $has_results -eq 0 ]] && out+=(--results /data)
  [[ $has_host -eq 0 ]] && out+=(--host 0.0.0.0)
  printf '%s\n' "${out[@]}"
}

case "${1:-}" in
  build)
    exec "${compose[@]}" build cli
    ;;
  dashboard)
    shift
    ensure_data
    mapfile -t dash_args < <(dashboard_args "$@")
    exec "${compose[@]}" --profile dashboard run --rm --service-ports dashboard "${dash_args[@]}"
    ;;
  shell)
    ensure_data
    exec "${compose[@]}" --profile cli run --rm --entrypoint bash cli
    ;;
  fetch)
    shift
    ensure_data
    exec "${compose[@]}" --profile cli run --rm --entrypoint python cli scripts/fetch_data.py --auto "$@"
    ;;
  --entrypoint)
    ep="$2"
    shift 2
    ensure_data
    exec "${compose[@]}" --profile cli run --rm --entrypoint "$ep" cli "$@"
    ;;
  *)
    ensure_data
    exec "${compose[@]}" --profile cli run --rm cli "$@"
    ;;
esac
