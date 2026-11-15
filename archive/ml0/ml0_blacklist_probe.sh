set -euo pipefail
ls -lh /data/ref/raw/ 2>/dev/null || true
ls -lh /data/ref/*.bed 2>/dev/null || true
# Try alternate fetch for blacklist with browser-like UA
python - <<'PY'
import urllib.request
from pathlib import Path
urls = [
  "https://cf.10xgenomics.com/supp/genome/GRCh38/svblacklist.bed",
  "https://github.com/10XGenomics/longranger/raw/master/longranger-cs/*/tenkit/reference/svblacklist.bed",
]
# Known copy often mirrored; try ENCODE hg38 blacklist as last resort only if 10x fails permanently
ua = "Mozilla/5.0 (compatible; arcsv2-fetch/1)"
for url in [
  "https://cf.10xgenomics.com/supp/genome/GRCh38/svblacklist.bed",
]:
  req = urllib.request.Request(url, headers={"User-Agent": ua, "Accept": "*/*"})
  try:
    with urllib.request.urlopen(req, timeout=60) as r:
      data = r.read(200)
      print("OK", url, "bytes_head", len(data), data[:80])
  except Exception as e:
    print("FAIL", url, type(e).__name__, e)
PY