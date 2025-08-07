# panget - a toolkit for Pangenome data processing

## Pangenome FASTA Merger with Header Annotation

This script merges multiple genome FASTA files into a single bgzipped FASTA file, while updating each sequence header to embed the sample name and haplotype ID. It also produces a per-genome scaffold mapping table that records how each sequence header was renamed.    
The script can be used to process fasta sequences for pangenomes.     

---

## What it does

- Accepts a list of FASTA files and their sample names.
- Renames all sequences using this format:  
```
SampleName#HaploID#OriginalHeader
```
- Concatenates all modified sequences into a single file.
- Compresses the output using `bgzip`.
- Writes a `.txt` file per input genome, listing the mapping between original and new scaffold names in tab-separated format.

---

## Requirements

Install the following Python and system dependencies:

```bash
./setup_pangenome_env.sh
```
That will install 
* Python biopython
* Python bgzip
* htslib for gzip

## Usage

```
usage: prepare-pangenomes.py [-h] [--haplo_id HAPLO_ID] [--threads THREADS]
                             genome_list output_fasta scaffold_map_dir

Concatenate input FASTAs into a bgzipped merged FASTA with modified headers.
Writes per-genome scaffold name maps.

positional arguments:
  genome_list          Input list: SAMPLE_NAME path/to/fasta per line
  output_fasta         Final output file (must end with .fasta.gz)
  scaffold_map_dir     Directory to save scaffold name mapping .txt files

optional arguments:
  -h, --help           show this help message and exit
  --haplo_id HAPLO_ID  Haplotype ID to embed in headers (default: 1)
  --threads THREADS    Threads to use for bgzip (default: 4)
```

Example
```
python prepare-pangenomes.py genome_list.txt \
	output.fasta.gz \
	scaffold_maps_dir\
	--haplo_id 1 \
	--threads 4
```

Positional arguments
| Argument             | Description                                               |
| -------------------- | --------------------------------------------------------- |
| `genome_list.txt`    | Path to input list file (see format below)                |
| `output.fasta.gz`    | Name for the final merged and bgzipped FASTA file         |
| `scaffold_maps_dir/` | Directory where per-genome `.txt` scaffold maps are saved |

Optional arguments
| Flag         | Description                                      | Default |
| ------------ | ------------------------------------------------ | ------- |
| `--haplo_id` | Haplotype ID to include in the header            | `1`     |
| `--threads`  | Number of threads to use for `bgzip` compression | `4`     |

## Input format

The input file is a tab- or space-delimited plain text file listing: 
```
SampleName    /path/to/sample1.fasta
SampleName2   /path/to/sample2.fasta.gz
```
FASTA files can be plain text or gzipped (.fasta or .fasta.gz).

## Output
* Merged FASTA
Final bgzipped file: `output.fasta.gz`

* Scaffold maps
One `.txt` file per genome in the given directory, named after the input file:
```
original_scaffold_id<TAB>new_scaffold_id
```
Example
```
scaffold1    Sample1#1#scaffold1
contig23     Sample1#1#contig23
```

## Example output structure
```
output.fasta.gz
scaffold_maps/
  ├── Sample1.txt
  ├── Sample2.txt

```

## Note on Haplotype Handling
This script is currently hardcoded to process only one haplotype per genome. The haplotype ID (--haplo_id) is applied uniformly to all sequences from each input FASTA file, and multiple haplotypes per sample are not yet supported. If you have multiple haplotypes per sample (e.g., phased assemblies), you'll need to run the script separately for each or modify the script to handle more complex input structures. 
Future versions may support multi-haplotype merging via an extended input format or optional flags.

## Licence 
MIT License

## Author and researcher
Kristina Gagalova      

Please get in touch with me if you have any issues or questions in running the script.   
