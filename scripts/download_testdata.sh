#!/usr/bin/env bash
# Download nf-core/test-datasets riboseq GSE182201 chr20 subset FASTQs into test/GSE182201/
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
INFO="${RIBO_SEQ_INFO_CSV:-$ROOT/test/info.csv}"
OUT="${1:-$ROOT/test/GSE182201}"
BASE="${RIBO_SEQ_TESTDATA_BASE:-https://raw.githubusercontent.com/nf-core/test-datasets/riboseq/testdata/GSE182201}"

if [ ! -f "$INFO" ]; then
  echo "Missing info.csv: $INFO" >&2
  exit 1
fi

mkdir -p "$OUT"
fetch() {
  local url="$1" dest="$2"
  if command -v wget >/dev/null 2>&1; then
    wget -c -O "$dest" "$url"
  else
    curl -L --fail -o "$dest" "$url"
  fi
}

awk -F, 'NR>1{
  for (i=2; i<=3; i++) {
    if ($i != "") print $i
  }
}' "$INFO" \
  | sort -u \
  | while read -r f; do
      [ -z "$f" ] && continue
      name="$(basename "$f")"
      dest="$OUT/$name"
      if [ -f "$dest" ]; then
        echo "skip existing $name"
        continue
      fi
      echo "fetch $name"
      fetch "$BASE/$name" "$dest"
    done

echo "Done. FASTQs in $OUT"
