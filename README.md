# AT_GC_compute

This code computes **GC and AT content maps** for a genome.  
It runs in a Python environment and requires:  

- A FASTA file (`genome.fa`)  
- A chromosome size file (`genome.chrom.sizes`)  

---

## Usage

Run the script as follows:

```bash
python compute_gc_at.py genome.fa genome.chrom.sizes 1kb genome_basecomp
