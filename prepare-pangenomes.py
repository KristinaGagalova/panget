#!/usr/bin/env python3

import argparse
import gzip
import sys
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

def load_checkpoint(path: Path):
    if not path.exists():
        return set()
    with open(path, "r") as f:
        done = {line.strip() for line in f if line.strip()}
    return done

def append_checkpoint(path: Path, sample_name: str):
    with open(path, "a") as f:
        f.write(sample_name + "\n")

def log_error(errlog: Path, sample_name: str, fasta_path: Path, exc: Exception):
    errlog.parent.mkdir(parents=True, exist_ok=True)
    with open(errlog, "a") as f:
        f.write(f"[ERROR] sample={sample_name} file={fasta_path} error={repr(exc)}\n")

def open_fasta_auto(p: Path):
    if p.suffix == ".gz":
        return gzip.open(p, "rt")
    return open(p, "r")

def process_one_genome(sample_name, fasta_path, out_handle, scaffold_map_dir: Path, delim="#", haplo_id="1"):
    """Process a single genome: rewrite headers -> write to out_handle; create scaffold map."""
    fasta_path = Path(fasta_path)
    # Write scaffold map file name (use base without .fasta if present)
    base_name = fasta_path.stem
    if base_name.endswith(".fasta"):  # handle .fasta.gz
        base_name = base_name[:-6]
    scaffold_map_dir.mkdir(parents=True, exist_ok=True)
    scaffold_map_path = scaffold_map_dir / f"{base_name}.txt"

    # Write to a temp map first; rename on success for atomicity
    tmp_map = scaffold_map_path.with_suffix(scaffold_map_path.suffix + ".tmp")

    with open(tmp_map, "w") as map_handle:
        with open_fasta_auto(fasta_path) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                new_id = f"{sample_name}{delim}{haplo_id}{delim}{record.id}"
                map_handle.write(f"{record.id}\t{new_id}\n")
                record.id = new_id
                record.description = ""
                SeqIO.write(record, out_handle, "fasta")

    tmp_map.rename(scaffold_map_path)  # atomic replace
    return scaffold_map_path

def rewrite_headers_and_concat_with_checkpoints(
    genomes, temp_output: Path, scaffold_map_dir: Path, ckpt_path: Path,
    errlog_path: Path, delim="#", haplo_id="1"
):
    """Rewrite headers + concatenate, skipping items in checkpoint; log errors and continue."""
    done = load_checkpoint(ckpt_path)
    total = len(genomes)
    to_do = [(s, p) for (s, p) in genomes if s not in done]

    if not temp_output.exists():
        # Start new file
        mode = "w"
    else:
        # Append if resuming
        mode = "a"

    print(f"[INFO] Genomes total: {total} | Completed (checkpoint): {len(done)} | Remaining: {len(to_do)}")
    failures = 0

    with open(temp_output, mode) as out_handle:
        for i, (sample_name, fasta_path) in enumerate(to_do, 1):
            print(f"[READ {i}/{len(to_do)}] {sample_name}: {fasta_path}")
            try:
                map_path = process_one_genome(sample_name, fasta_path, out_handle, scaffold_map_dir, delim, haplo_id)
                append_checkpoint(ckpt_path, sample_name)
                print(f"[WRITE] Scaffold map saved: {map_path} | Checkpointed: {sample_name}")
            except Exception as e:
                failures += 1
                log_error(errlog_path, sample_name, Path(fasta_path), e)
                print(f"[WARN] Failed {sample_name}; logged and continuing.")

    print(f"[WRITE] Temporary FASTA written/appended: {temp_output}")
    return failures

def bgzip_file(input_path, threads=4):
    """Compress the file using bgzip"""
    print(f"[BGZIP] Compressing {input_path} with {threads} threads...")
    subprocess.run(["bgzip", "-f", "-@", str(threads), str(input_path)], check=True)
    print(f"[DONE] Compressed to: {input_path}.gz")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Concatenate input FASTAs into a bgzipped merged FASTA with modified headers. Writes per-genome scaffold name maps. Supports checkpointed resume."
    )
    parser.add_argument("genome_list", help="Input list: SAMPLE_NAME path/to/fasta per line")
    parser.add_argument("output_fasta", help="Final output file (must end with .fasta.gz)")
    parser.add_argument("scaffold_map_dir", help="Directory to save scaffold name mapping .txt files")
    parser.add_argument("--haplo_id", default="1", help="Haplotype ID to embed in headers (default: 1)")
    parser.add_argument("--threads", type=int, default=4, help="Threads to use for bgzip (default: 4)")
    parser.add_argument("--checkpoint", type=Path, help="Path to checkpoint file (defaults to <output>.ckpt)")
    parser.add_argument("--error_log", type=Path, help="Path to error log file (defaults to <output>.errors.log)")
    parser.add_argument("--reset", action="store_true", help="Ignore existing checkpoint and start fresh (overwrites temp output)")
    args = parser.parse_args()

    output_path = Path(args.output_fasta)
    if output_path.suffix != ".gz":
        raise ValueError("Output file must end with '.gz' for bgzip compression.")

    temp_path = output_path.with_suffix('')  # uncompressed temp (same name without .gz)
    ckpt_path = args.checkpoint or (output_path.parent / (output_path.name + ".ckpt"))
    errlog_path = args.error_log or (output_path.parent / (output_path.name + ".errors.log"))

    if args.reset:
        if ckpt_path.exists():
            ckpt_path.unlink()
        if temp_path.exists():
            temp_path.unlink()
        print("[INFO] Reset requested: removed existing checkpoint and temp output.")

    genome_entries = read_genome_list(args.genome_list)
    failures = rewrite_headers_and_concat_with_checkpoints(
        genome_entries,
        temp_output=temp_path,
        scaffold_map_dir=Path(args.scaffold_map_dir),
        ckpt_path=ckpt_path,
        errlog_path=errlog_path,
        haplo_id=args.haplo_id
    )

    # Compress even if there were failures—your pipeline may still want partial results.
    bgzip_file(temp_path, threads=args.threads)

    if failures > 0:
        print(f"[DONE] Completed with {failures} failures. See log: {errlog_path}")
        sys.exit(1)
    else:
        print("[DONE] All genomes processed successfully.")
