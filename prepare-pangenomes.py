#!/usr/bin/env python3

import argparse
import gzip
import subprocess
from pathlib import Path
from Bio import SeqIO

def read_genome_list(file_path):
    """Read input list with: SAMPLE_NAME path/to/file.fasta"""
    genomes = []
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split(None, 1)
            if len(parts) == 2:
                name, fasta_path = parts
                genomes.append((name, fasta_path))
    return genomes

def rewrite_headers_and_concat(genomes, temp_output, scaffold_map_dir, delim="#", haplo_id="1"):
    """Rewrite headers and concatenate all FASTA files into a temporary output (handles .gz input)"""
    genome_count = len(genomes)
    print(f"[INFO] Preparing to add {genome_count} genomes to the pangenome")

    scaffold_map_dir.mkdir(parents=True, exist_ok=True)

    with open(temp_output, "w") as out_handle:
        for sample_name, fasta_path in genomes:
            fasta_path = Path(fasta_path)
            print(f"[READ] {fasta_path}")

            # Choose open function based on file extension
            if fasta_path.suffix == ".gz":
                open_func = gzip.open
                mode = "rt"
            else:
                open_func = open
                mode = "r"

            # Write scaffold map file
            base_name = fasta_path.stem
            if base_name.endswith(".fasta"):  # handle .fasta.gz
                base_name = base_name[:-6]
            scaffold_map_path = scaffold_map_dir / f"{base_name}.txt"

            with open(scaffold_map_path, "w") as map_handle:
                with open_func(fasta_path, mode) as handle:
                    for record in SeqIO.parse(handle, "fasta"):
                        new_id = f"{sample_name}{delim}{haplo_id}{delim}{record.id}"
                        map_handle.write(f"{record.id}\t{new_id}\n")

                        record.id = new_id
                        record.description = ""
                        SeqIO.write(record, out_handle, "fasta")

            print(f"[WRITE] Scaffold map saved: {scaffold_map_path}")

    print(f"[WRITE] Temporary FASTA written: {temp_output}")

def bgzip_file(input_path, threads=4):
    """Compress the file using bgzip"""
    print(f"[BGZIP] Compressing {input_path} with {threads} threads...")
    subprocess.run(["bgzip", "-f", "-@", str(threads), str(input_path)], check=True)
    print(f"[DONE] Compressed to: {input_path}.gz")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Concatenate input FASTAs into a bgzipped merged FASTA with modified headers. Writes per-genome scaffold name maps."
    )
    parser.add_argument("genome_list", help="Input list: SAMPLE_NAME path/to/fasta per line")
    parser.add_argument("output_fasta", help="Final output file (must end with .fasta.gz)")
    parser.add_argument("scaffold_map_dir", help="Directory to save scaffold name mapping .txt files")
    parser.add_argument("--haplo_id", default="1", help="Haplotype ID to embed in headers (default: 1)")
    parser.add_argument("--threads", type=int, default=4, help="Threads to use for bgzip (default: 4)")
    args = parser.parse_args()

    output_path = Path(args.output_fasta)
    if output_path.suffix != ".gz":
        raise ValueError("Output file must end with '.gz' for bgzip compression.")

    temp_path = output_path.with_suffix('')  # uncompressed temp

    genome_entries = read_genome_list(args.genome_list)
    rewrite_headers_and_concat(
        genome_entries,
        temp_output=temp_path,
        scaffold_map_dir=Path(args.scaffold_map_dir),
        haplo_id=args.haplo_id
    )
    bgzip_file(temp_path, threads=args.threads)
