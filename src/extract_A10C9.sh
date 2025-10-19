#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./extract_A10C09.sh [--rename] <genome.(fa|fasta)[.gz]> [output.fa]
# Examples:
#   ./extract_A10C09.sh genome.fasta
#   ./extract_A10C09.sh --rename genome.fa.gz
#   ./extract_A10C09.sh --rename genome.fasta subset.fa

module load samtools/1.15--h3843a85_0 2>/dev/null || true

SUFFIX="A10C09"
RENAME=0

# ---- parse options ----
args=()
while (( "$#" )); do
  case "$1" in
    --rename) RENAME=1; shift ;;
    -h|--help)
      echo "Usage: $0 [--rename] <genome.(fa|fasta)[.gz]> [output.fa]"
      exit 0
      ;;
    --) shift; break ;;
    -*) echo "Unknown option: $1" >&2; exit 1 ;;
    *)  args+=("$1"); shift ;;
  esac
done
set -- "${args[@]}"

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 [--rename] <genome.(fa|fasta)[.gz]> [output.fa]" >&2
  exit 1
fi

FASTA="$1"

# ---- robust base name extraction (handles .fa/.fasta with optional .gz) ----
# Strip one .gz if present, then strip .fa or .fasta (in that order)
NO_GZ="${FASTA%.gz}"
BASE="$(basename "${NO_GZ%.fa}")"
BASE="$(basename "${BASE%.fasta}")"

# If user provided output, use it; else force .fa extension with single suffix
if [[ $# -eq 2 ]]; then
  OUT="$2"
else
  OUT="${BASE}.${SUFFIX}.fa"
fi

# Skip if input already looks like a subset to avoid reprocessing
# (match either SUFFIX or the older A10C9 just in case)
if [[ "$BASE" =~ \.${SUFFIX}$ || "$BASE" =~ \.A10C9$ ]]; then
  echo "⚠️  Input already appears to be a subset ($BASE). Skipping $FASTA."
  exit 0
fi

# Early exit if output exists
if [[ -f "$OUT" ]]; then
  echo "⚠️  Output exists: $OUT — skipping $FASTA."
  exit 0
fi

FAI="${FASTA}.fai"
if [[ ! -f "$FAI" ]]; then
  echo "Indexing $FASTA ..."
  samtools faidx "$FASTA"
fi

# Build list of target IDs and a name→token map
TMP_IDS="$(mktemp)"
TMP_MAP="$(mktemp)"   # orig_token \t TOKEN(A/Cxx)

awk '
  BEGIN{ IGNORECASE=1 }
  /^>/ {
    h=$1; sub(/^>/,"",h)                 # first token (samtools uses up to whitespace)
    if (match(h, /(A0[1-9]|A10|C0[1-9])/)) {
      token=substr(h, RSTART, RLENGTH)
      print h > "'"$TMP_IDS"'"
      print h "\t" token > "'"$TMP_MAP"'"
    }
  }
' "$FASTA"

sort -u -o "$TMP_IDS" "$TMP_IDS"
sort -u -o "$TMP_MAP" "$TMP_MAP"

if [[ ! -s "$TMP_IDS" ]]; then
  echo "No matching scaffolds (A01–A10 or C01–C09) in $FASTA."
  rm -f "$TMP_IDS" "$TMP_MAP"
  exit 2
fi

echo "Extracting $(wc -l < "$TMP_IDS") scaffolds to $OUT ..."
samtools faidx "$FASTA" -r "$TMP_IDS" > "$OUT"

if [[ "$RENAME" -eq 1 ]]; then
  echo "Renaming headers (Axx/Cxx) ..."
  awk -v mapfile="$TMP_MAP" '
    BEGIN { FS=OFS="\t"; while((getline<mapfile)>0) map[$1]=$2; close(mapfile) }
    /^>/ {
      name=$1; sub(/^>/,"",name)
      if (name in map) print ">" map[name]; else print $0
      next
    }
    { print }
  ' "$OUT" > "${OUT}.renamed.tmp"
  mv "${OUT}.renamed.tmp" "$OUT"
fi

rm -f "$TMP_IDS" "$TMP_MAP"
echo "✅ Done: $OUT"
