from pathlib import Path
import hashlib

ROOT = Path(__file__).resolve().parents[1]
manifest = ROOT / "SHA256SUMS.txt"

bad = []
for line in manifest.read_text(encoding="utf-8").splitlines():
    if not line.strip():
        continue
    expected, rel = line.split("  ", 1)
    p = ROOT / rel
    if not p.exists():
        bad.append((rel, "missing"))
        continue
    h = hashlib.sha256(p.read_bytes()).hexdigest()
    if h != expected:
        bad.append((rel, h))

if bad:
    for item in bad:
        print("FAIL", *item)
    raise SystemExit(1)
print("All release checksums verified.")
