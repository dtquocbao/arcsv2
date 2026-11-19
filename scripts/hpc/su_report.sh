#!/usr/bin/env bash
# Usage: bash scripts/hpc/su_report.sh <snakemake log>
#
# Reads SLURM job ids from a Snakemake log (one "submitted with SLURM jobid"
# line per job) and asks sacct for
#   JobID,JobName,Elapsed,AllocCPUS,ReqMem,MaxRSS
# Service units use the Grace rule: elapsed hours times
# max(AllocCPUS, requested GB / 7.5). ReqMem suffixes: M = mebibytes, G = gibibytes,
# and a trailing c means per-CPU. "SU per genome" is fetch_region + call + merge.
set -euo pipefail

if [[ $# -ne 1 || ! -f "${1:-}" ]]; then
  echo "usage: $0 <snakemake log>" >&2
  exit 2
fi
LOG="$1"

declare -A SMK_RULE=()
declare -A RULE_IDS=()
current_rule=""

while IFS= read -r line || [[ -n "$line" ]]; do
  if [[ "$line" =~ ^rule[[:space:]]+([^:]+): ]]; then
    current_rule="${BASH_REMATCH[1]}"
    continue
  fi
  if [[ "$line" =~ ^[[:space:]]*jobid:[[:space:]]*([0-9]+) ]]; then
    if [[ -n "$current_rule" ]]; then
      SMK_RULE["${BASH_REMATCH[1]}"]="$current_rule"
    fi
    continue
  fi
  smk_id=""
  slurm_id=""
  if [[ "$line" =~ Job[[:space:]]+([0-9]+)[[:space:]]+has[[:space:]]+been[[:space:]]+submitted[[:space:]]+with[[:space:]]+SLURM[[:space:]]+jobid[[:space:]]+([0-9]+) ]]; then
    smk_id="${BASH_REMATCH[1]}"
    slurm_id="${BASH_REMATCH[2]}"
  elif [[ "$line" =~ external[[:space:]]+jobid[[:space:]]+[\'\"]?([0-9]+) ]]; then
    slurm_id="${BASH_REMATCH[1]}"
  fi
  if [[ -n "$slurm_id" ]]; then
    rule="${current_rule:-unknown}"
    if [[ -n "$smk_id" && -n "${SMK_RULE[$smk_id]:-}" ]]; then
      rule="${SMK_RULE[$smk_id]}"
    fi
    RULE_IDS["$rule"]+="${slurm_id} "
  fi
done < "$LOG"

if [[ ${#RULE_IDS[@]} -eq 0 ]]; then
  echo "No SLURM job ids in ${LOG}" >&2
  exit 1
fi

ids=()
for rule in "${!RULE_IDS[@]}"; do
  for id in ${RULE_IDS[$rule]}; do
    ids+=("$id")
  done
done

if ! command -v sacct >/dev/null 2>&1; then
  echo "sacct is not on PATH. Run this on a Grace login node." >&2
  exit 1
fi

# -P is the parsable form of the requested columns. -n drops the header.
# Parent rows only: step ids (123.batch) would double-count the allocation.
sacct_out="$(sacct -j "$(IFS=,; echo "${ids[*]}")" \
  --format=JobID,JobName,Elapsed,AllocCPUS,ReqMem,MaxRSS -P -n)"

elapsed_hours() {
  local t="$1" d=0 h=0 m=0 s=0 a b c
  if [[ "$t" == *-* ]]; then
    d="${t%%-*}"
    t="${t#*-}"
  fi
  IFS=: read -r a b c <<< "$t"
  if [[ -z "${c:-}" ]]; then
    m="${a:-0}"
    s="${b:-0}"
  else
    h="${a:-0}"
    m="${b:-0}"
    s="${c:-0}"
  fi
  awk -v d="$d" -v h="$h" -v m="$m" -v s="$s" 'BEGIN { printf "%.6f", (d * 24) + h + (m / 60) + (s / 3600) }'
}

# Print one SU for a sacct ReqMem token and an AllocCPUS count.
su_for() {
  local elapsed="$1" cpus="$2" req="$3"
  local hours num unit per_cpu=0 gb
  hours="$(elapsed_hours "$elapsed")"
  num="$(printf '%s' "$req" | sed -E 's/^([0-9.]+).*/\1/')"
  unit="$(printf '%s' "$req" | sed -E 's/^[0-9.]+//; s/[+]//g')"
  if [[ "$unit" == *c ]]; then
    per_cpu=1
    unit="${unit%c}"
  fi
  case "$unit" in
    K|k) gb="$(awk -v n="$num" 'BEGIN { printf "%.6f", n / 1024 / 1024 }')" ;;
    M|m|"") gb="$(awk -v n="${num:-0}" 'BEGIN { printf "%.6f", n / 1024 }')" ;;
    G|g) gb="$(awk -v n="$num" 'BEGIN { printf "%.6f", n }')" ;;
    T|t) gb="$(awk -v n="$num" 'BEGIN { printf "%.6f", n * 1024 }')" ;;
    *) gb="$(awk -v n="${num:-0}" 'BEGIN { printf "%.6f", n / 1024 }')" ;;
  esac
  if [[ "$per_cpu" -eq 1 ]]; then
    gb="$(awk -v g="$gb" -v c="$cpus" 'BEGIN { printf "%.6f", g * c }')"
  fi
  awk -v h="$hours" -v c="${cpus:-0}" -v g="$gb" 'BEGIN {
    mem = g / 7.5
    charge = c + 0
    if (mem > charge) charge = mem
    printf "%.4f", h * charge
  }'
}

declare -A RULE_SU=()
declare -A RULE_N=()
echo "rule jobs SU elapsed_example"
while IFS='|' read -r jobid _jobname elapsed cpus reqmem _maxrss; do
  [[ -z "${jobid:-}" ]] && continue
  [[ "$jobid" == *.* ]] && continue
  rule=""
  for name in "${!RULE_IDS[@]}"; do
    for id in ${RULE_IDS[$name]}; do
      if [[ "$id" == "$jobid" ]]; then
        rule="$name"
      fi
    done
  done
  [[ -z "$rule" ]] && continue
  su="$(su_for "$elapsed" "${cpus:-0}" "${reqmem:-0}")"
  RULE_SU["$rule"]="$(awk -v a="${RULE_SU[$rule]:-0}" -v b="$su" 'BEGIN { printf "%.4f", a + b }')"
  RULE_N["$rule"]=$(( ${RULE_N[$rule]:-0} + 1 ))
done <<< "$sacct_out"

genome="0"
for rule in $(printf '%s\n' "${!RULE_SU[@]}" | sort); do
  printf '%s %s %s\n' "$rule" "${RULE_N[$rule]}" "${RULE_SU[$rule]}"
  if [[ "$rule" == "fetch_region" || "$rule" == "call" || "$rule" == "merge" ]]; then
    genome="$(awk -v a="$genome" -v b="${RULE_SU[$rule]}" 'BEGIN { printf "%.4f", a + b }')"
  fi
done
printf 'SU per genome %s\n' "$genome"
