#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./rename_fasta.sh <map.txt> <genome.(fa|fasta|fna)[.gz]>
# Output:
#   genome.renamed.fa
# Map format:
#   First token  = old ID (e.g. CP123876.1)
#   Last  token  = new name (e.g. A01)
#   Middle text ignored; may contain spaces/tabs.
#   Example:
#   CP123876.1 chromosome[TAB]A01
#   .....

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <map.txt> <genome.(fa|fasta|fna)[.gz]>" >&2
  exit 1
fi

MAP="$1"
FA="$2"

[[ -f "$MAP" ]] || { echo "❌ Mapping file not found: $MAP" >&2; exit 2; }
[[ -f "$FA"  ]] || { echo "❌ FASTA file not found: $FA" >&2; exit 3; }

# Output name
NO_GZ="${FA%.gz}"
BASE="$(basename "${NO_GZ%.*}")"
OUT="${BASE}.renamed.fa"

echo "Renaming headers in $FA using $MAP → $OUT"

# Handle gzip automatically
if [[ "$FA" == *.gz ]]; then
  DECOMPRESS="gzip -dc -- \"$FA\""
else
  DECOMPRESS="cat -- \"$FA\""
fi

# Run AWK
eval "$DECOMPRESS" | awk -v mapfile="$MAP" '
  BEGIN {
    FS="[ \t]+"; OFS="\t"
    # load map: first token -> last token
    while ((getline line < mapfile) > 0) {
      gsub(/^[ \t]+|[ \t]+$/, "", line)
      if (line=="" || line ~ /^#/) continue
      n = split(line, f, /[ \t]+/)
      if (n >= 2) {
        oldid = f[1]; newid = f[n]
        m[oldid] = newid
      }
    }
    close(mapfile)
  }
  /^>/ {
    hdr=$0; sub(/^>/,"",hdr)
    split(hdr,a,/[\t ]+/); id=a[1]
    if (id in m) print ">" m[id]; else print $0
    next
  }
  { print }
' > "$OUT"

echo "✅ Done: $OUT"
