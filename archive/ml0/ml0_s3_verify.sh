set -euo pipefail
export PYTHONPATH=
for s in HG01884 HG01150 HG00126; do
  samtools quickcheck /data/samples/$s/reads/$s.cram && echo "quickcheck_ok $s"
  samtools idxstats /data/samples/$s/reads/$s.cram | awk -v s=$s '$1=="chr20"{print s,$1,$2,$3}'
done
test -f /data/ref/simple_repeats.bed && test -f /data/ref/segdups.bed && echo filter_beds_ok
snakemake -s workflows/label.smk --configfile workflows/config/label.ml0_s3.yaml -n --quiet 2>&1 | tail -20