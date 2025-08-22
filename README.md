# AT_GC_compute
The code can be used to compute the GC and AT map for a genome
This code can be executed in the python environment and requires a fasta file and chromosome size file:
"python compute_gc_at.py genome.fa genome.chrom.sizes 1kb genome_basecomp"
Here, compute_gc_at.py is the code part, genome.fa could be any genome, genome.chrom.size and 1kb that is window size (it can be changed as we wish), and the output name.
To compute the genome chr.size, it can be done as following:
Make a fai index "samtools faidx genome.fa"
then use cut command "cut -f1,2 genome.fa.fai > genome.chrom.sizes"
