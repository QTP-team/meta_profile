import argparse
import os
from Bio import SeqIO

parser = argparse.ArgumentParser(
        usage = 'python genome_len.py -i *.fa -o genome_contig_length.txt -s .fa')
parser.add_argument("-i", "--infiles", nargs = "+", help = "One or more genomes")
parser.add_argument("-o", "--output")
parser.add_argument("-s", "--suffix", default = ".fa", help = "suffix of genome files")

args = parser.parse_args()
f_o = open(args.output, "w")
f_o.write("Genome_ID\tContig_ID\tContig_length\n")

n_genome = 0
n_contig = 0
for fa in args.infiles:
    gID = os.path.basename(fa).replace(args.suffix, "")
    n_genome += 1
    for seq_r in SeqIO.parse(fa, "fasta"):
        f_o.write("%s\t%s\t%s\n" %(gID, seq_r.id, len(seq_r)))
        n_contig += 1

print("Number of genome: %s, Number of contig: %s" %(n_genome, n_contig))
