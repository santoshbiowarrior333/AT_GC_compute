#!/usr/bin/env python3
import sys
from pyfaidx import Fasta

def parse_window_size(window_str):
    """Convert window size input (e.g. 1kb, 50kb, 1000) into integer bp"""
    window_str = window_str.lower()
    if window_str.endswith("kb"):
        return int(float(window_str[:-2]) * 1000)
    elif window_str.endswith("mb"):
        return int(float(window_str[:-2]) * 1_000_000)
    else:
        return int(window_str)

def compute_gc_at(fasta_file, chrom_sizes_file, window_size, out_prefix):
    # Load genome
    print(f"Loading genome from {fasta_file} ...")
    genome = Fasta(fasta_file, sequence_always_upper=True)
    print(f"Genome loaded with {len(genome)} chromosomes")

    # Load chromosome sizes
    chrom_sizes = {}
    with open(chrom_sizes_file) as f:
        for line in f:
            chrom, size = line.strip().split()
            chrom_sizes[chrom] = int(size)
    print(f"Loaded chromosome sizes for {len(chrom_sizes)} chromosomes")

    gc_out = open(f"{out_prefix}_GC.bedGraph", "w")
    at_out = open(f"{out_prefix}_AT.bedGraph", "w")

    for chrom, size in chrom_sizes.items():
        if chrom not in genome:
            print(f"Warning: {chrom} not found in FASTA, skipping")
            continue
        print(f"Processing chromosome {chrom}, size {size} bp")
        seq = genome[chrom][:].seq

        num_windows = (size + window_size - 1) // window_size
        for i, start in enumerate(range(0, size, window_size), 1):
            end = min(start + window_size, size)
            window_seq = seq[start:end]

            a = window_seq.count("A")
            t = window_seq.count("T")
            g = window_seq.count("G")
            c = window_seq.count("C")

            total = a + t + g + c
            if total == 0:
                continue

            at = (a + t) / total
            gc = (g + c) / total

            gc_out.write(f"{chrom}\t{start}\t{end}\t{gc:.4f}\n")
            at_out.write(f"{chrom}\t{start}\t{end}\t{at:.4f}\n")

            # Print progress every 1000 windows or at the last window
            if i % 1000 == 0 or i == num_windows:
                print(f"  Window {i}/{num_windows} ({start}-{end} bp) completed")

    gc_out.close()
    at_out.close()
    print(f"Finished. Output files: {out_prefix}_GC.bedGraph, {out_prefix}_AT.bedGraph")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: compute_gc_at.py <reference.fasta> <chrom.sizes> <window_size> <out_prefix>")
        print("Example: compute_gc_at.py genome.fa genome.chrom.sizes 1kb genome_basecomp")
        sys.exit(1)

    fasta_file = sys.argv[1]
    chrom_sizes_file = sys.argv[2]
    window_size = parse_window_size(sys.argv[3])
    out_prefix = sys.argv[4]

    compute_gc_at(fasta_file, chrom_sizes_file, window_size, out_prefix)
