#!/usr/bin/env bash
set -euo pipefail

module load parallel/20240822
module load samtools/1.15--h3843a85_0

INDEX=/path/to/index/all_genomes.fasta.gz.fai
GENOMES=/path/to/genomes/all_genomes.fasta.gz
NAMES=names.txt   # one name/pattern per line (e.g., A10, chrC03, etc.)
# NAMES contains a list of chromosome names to be parsed

export INDEX GENOMES

# sanity: ensure names are non-empty, ignore comment lines
grep -v -E '^\s*(#|$)' "$NAMES" | parallel --env GENOMES --env INDEX -j 12 --linebuffer '
  NAME="{}"
  OUT_MATCH="group.${NAME}.fasta.gz"
  OUT_REST="group.${NAME}.rest.fasta.gz"

  MATCH_IDS=$(mktemp)
  REST_IDS=$(mktemp)
  ALL_IDS=$(mktemp)

  # Collect all IDs (column 1 of .fai)
  awk -F"\t" '\''{print $1}'\'' "$INDEX" | sort -u > "$ALL_IDS"

  # Fixed substring match on the 3rd "#" field (case-sensitive).
  # To make it case-insensitive, wrap both sides with tolower().
  awk -F"\t" -v n="$NAME" '\''{
      split($1,a,"#");
      if (index(a[3], n) > 0) print $1
  }'\'' "$INDEX" | sort -u > "$MATCH_IDS"

  # Everything else = ALL minus MATCH
  comm -23 "$ALL_IDS" "$MATCH_IDS" > "$REST_IDS"

  if [ -s "$MATCH_IDS" ]; then
    samtools faidx "$GENOMES" -r "$MATCH_IDS" | bgzip -@ 4 > "$OUT_MATCH"
    echo "Generated $OUT_MATCH"
  else
    echo "No matches for pattern: ${NAME}" >&2
  fi

  if [ -s "$REST_IDS" ]; then
    samtools faidx "$GENOMES" -r "$REST_IDS" | bgzip -@ 4 > "$OUT_REST"
    echo "Generated $OUT_REST"
  else
    echo "No REST sequences for pattern: ${NAME}" >&2
  fi

  rm -f "$MATCH_IDS" "$REST_IDS" "$ALL_IDS"
'
