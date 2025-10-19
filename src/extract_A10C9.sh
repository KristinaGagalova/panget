#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./extract_A10C09.sh [--rename] <genome.(fa|fasta|fna)[.gz]> [output.fa]
# Default output: NAME.A10C09.fa
# Behavior:
#   - Matches A01–A10 / C01–C09 and A1–A10 / C1–C9 anywhere in the header.
#   - By default, keeps FULL original headers in the output.
#   - With --rename, headers become >A01..A10 / >C01..C09.

module load samtools/1.15--h3843a85_0 2>/dev/null || true

SUFFIX="A10C9"
RENAME=0

# ---- parse options ----
args=()
while (( "$#" )); do
  case "$1" in
    --rename) RENAME=1; shift ;;
    -h|--help)
      echo "Usage: $0 [--rename] <genome.(fa|fasta|fna)[.gz]> [output.fa]"
      exit 0 ;;
    --) shift; break ;;
    -*) echo "Unknown option: $1" >&2; exit 1 ;;
    *)  args+=("$1"); shift ;;
  esac
done
set -- "${args[@]}"

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 [--rename] <genome.(fa|fasta|fna)[.gz]> [output.fa]" >&2
  exit 1
fi

FASTA="$1"

# ---- base name extraction (fa/fasta/fna + optional .gz) ----
NO_GZ="${FASTA%.gz}"
BASE="$(basename "${NO_GZ%.fa}")"
BASE="$(basename "${BASE%.fasta}")"
BASE="$(basename "${BASE%.fna}")"

# Output path
OUT="${2:-${BASE}.${SUFFIX}.fa}"

# Skip already-processed inputs / outputs
if [[ "$BASE" =~ \.${SUFFIX}$ || "$BASE" =~ \.A10C9$ ]]; then
  echo "⚠️  Input already appears to be a subset ($BASE). Skipping $FASTA."
  exit 0
fi
if [[ -f "$OUT" ]]; then
  echo "⚠️  Output exists: $OUT — skipping $FASTA."
  exit 0
fi

# Index if needed
FAI="${FASTA}.fai"
if [[ ! -f "$FAI" ]]; then
  echo "Indexing $FASTA ..."
  samtools faidx "$FASTA"
fi

# ---- Build ID lists and maps ----
TMP_IDS="$(mktemp)"
TMP_MAP="$(mktemp)"   # id \t normalized_token \t full_header_without_">"

awk '
  BEGIN { IGNORECASE=1 }
  /^>/ {
    full=$0; sub(/^>/,"",full)
    split(full, arr, /[ \t]/); id=arr[1]

    # Search entire header for A/C + (1..10) with optional leading zero
    if (match(full, /([aAcC])(0?[1-9]|10)/, m)) {
      letter=toupper(m[1]); num=m[2]+0
      if ((letter=="A" && num>=1 && num<=10) || (letter=="C" && num>=1 && num<=9)) {
        token=sprintf("%s%02d", letter, num)
        print id > "'"$TMP_IDS"'"
        print id "\t" token "\t" full > "'"$TMP_MAP"'"
      }
    }
  }
' "$FASTA"

sort -u -o "$TMP_IDS" "$TMP_IDS"
sort -u -o "$TMP_MAP" "$TMP_MAP"

if [[ ! -s "$TMP_IDS" ]]; then
  echo "No matching scaffolds (A1–A10 / C1–C9) in $FASTA."
  rm -f "$TMP_IDS" "$TMP_MAP"
  exit 2
fi

COUNT=$(wc -l < "$TMP_IDS")
echo "Extracting ${COUNT} scaffolds to $OUT ..."
samtools faidx "$FASTA" -r "$TMP_IDS" > "$OUT"

# ---- Rewrite headers: full original (default) or normalized tokens (with --rename) ----
# samtools faidx outputs headers as >ID only; we restore/rename here.
if [[ "$RENAME" -eq 1 ]]; then
  echo "Renaming headers to tokens (A01..A10 / C01..C09) ..."
  awk -v mapfile="$TMP_MAP" '
    BEGIN {
      FS=OFS="\t"
      while ((getline < mapfile) > 0) { id=$1; tok=$2; full[$1]=tok }
      close(mapfile)
    }
    /^>/ {
      id=$1; sub(/^>/,"",id)
      if (id in full) print ">" full[id]; else print $0
      next
    }
    { print }
  ' "$OUT" > "${OUT}.tmp" && mv "${OUT}.tmp" "$OUT"
else
  echo "Restoring full original headers ..."
  awk -v mapfile="$TMP_MAP" '
    BEGIN {
      FS=OFS="\t"
      while ((getline < mapfile) > 0) { id=$1; hdr=$3; full[id]=hdr }
      close(mapfile)
    }
    /^>/ {
      id=$1; sub(/^>/,"",id)
      if (id in full) print ">" full[id]; else print $0
      next
    }
    { print }
  ' "$OUT" > "${OUT}.tmp" && mv "${OUT}.tmp" "$OUT"
fi

rm -f "$TMP_IDS" "$TMP_MAP"
echo "✅ Done: $OUT"
