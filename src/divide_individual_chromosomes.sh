#!/usr/bin/bash

module load parallel/20240822
module load samtools/1.15--h3843a85_0

INDEX=/spath/to/index/all_genomes.fasta.gz.fai
GENOMES=/spath/to/genomes/all_genomes.fasta.gz

export INDEX GENOMES

awk -F'\t' '{split($1,a,"#"); if (a[3]!="") print a[3]}' "$INDEX" | sort -u |
parallel --env GENOMES --env INDEX -j 12 '
  CHROM="{}"
  CHR_FASTA="bnapus16.${CHROM}.fasta.gz"

  # exact match: take rows where 3rd "#" field == $CHROM, output full IDs (col1)
  samtools faidx "$GENOMES" -r <(
    awk -v c="$CHROM" -F"\t" '\''{split($1,a,"#"); if (a[3]==c) print $1}'\'' "$INDEX"
  ) | bgzip -@ 4 > "$CHR_FASTA"

  echo "Generated $CHR_FASTA"
'
