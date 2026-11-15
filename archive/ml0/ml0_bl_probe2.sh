set -euo pipefail
python - <<'PY'
import urllib.request
from pathlib import Path
candidates = [
  "https://raw.githubusercontent.com/10XGenomics/longranger/v2.2.2/longranger-cs/2.2.2/tenkit/lib/python/tenkit/reference.py",
  "https://github.com/10XGenomics/supernova/raw/master/tenkit/lib/python/tenkit/reference/sv_blacklist.bed",
  "https://raw.githubusercontent.com/igordot/genomics/master/data/10x_sv_blacklist_grch38.bed",
  "https://s3-us-west-2.amazonaws.com/10x.files/supp/genome/GRCh38/svblacklist.bed",
  "https://cdn.10xgenomics.com/supp/genome/GRCh38/svblacklist.bed",
]
ua = "Mozilla/5.0"
for url in candidates:
  req = urllib.request.Request(url, headers={"User-Agent": ua})
  try:
    with urllib.request.urlopen(req, timeout=30) as r:
      data = r.read(120)
      print("OK", r.status, url, data[:60])
  except Exception as e:
    print("FAIL", url, type(e).__name__, str(e)[:120])
PY