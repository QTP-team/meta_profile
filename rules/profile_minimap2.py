import argparse
import subprocess
import os

parser = argparse.ArgumentParser(
    prog = "profile_minimap2.py",
    usage = "python rules/profile_minimap2.py -x 1.index/merge_sp3_minimap2_index.min -r1 test1_data_R1.fastq -r2 test1_data_R2.fastq -contig 1.index/genome_contig_length.txt -t 8 -i 0.95 -n 100 -o test1")

arg = parser.add_argument

arg('-x', '--minimap2_index', type=str, help = 'The minimap2 index of the reference genome')
arg('-r1', '--rmhost_reads1', type=str, help = 'Query input files are FASTQ .fq/.fastq')
arg('-r2', '--rmhost_reads2', type=str, help = 'Query input files are FASTQ .fq/.fastq')
arg('-contig', '--genome_contig_length', type=str, help = 'The file of the reference genomes corresponding to the length of contigs')
arg('-t', '--threads', metavar = 'N', type=int, default=4, help = 'Threads for minimap2')
arg('-id', '--identity', metavar = 'N', type=float, default=0.95, help = 'Alignments which identity larger than N will be retained')
arg('-n', '--reads_num', metavar = 'N', type=int, default=100, help = 'Genomes which supported by at least N paired-end reads will be retained')
arg('-o', '--output_prefix', type=str, help = 'The prefix of the output result')


def profile(index,rmhost1,rmhost2,genome_contig_length_file,threads,identity,reads_num,output_prefix):
    os.mkdir(output_prefix)
    map_file = os.path.join(output_prefix, output_prefix + ".paf")
    paf_file = open(map_file, "w")
    cal_file = os.path.join(output_prefix, output_prefix + "_tmp.tsv")
    aligned_file = open(cal_file, "w")
    subprocess.run("minimap2 -cx sr -t %s %s %s %s -N 2 | awk '{print $1,$6,$8,$10,$11}'" %(threads,index,rmhost1,rmhost2),shell=True,encoding="utf-8",stdout=paf_file)
    paf_file.close()
    aligned_file.write("Read_ID\tContig_ID\tR1_start_loc\tR2_start_loc\tIdentity\n")
    with open(map_file, "r") as f:
        for line in f:
            r_ID, r_contig, r_start_loc, r_match_bases, r_alignment_block_length = line.strip().split(" ")
            if (r_ID.endswith("/1")):
                r1_ID, r1_contig, r1_start_loc, r1_match_bases, r1_alignment_block_length = r_ID, r_contig, r_start_loc, r_match_bases, r_alignment_block_length
            else:
                r2_ID, r2_contig, r2_start_loc, r2_match_bases, r2_alignment_block_length = r_ID, r_contig, r_start_loc, r_match_bases, r_alignment_block_length
                if r1_ID[:-2] == r2_ID[:-2] and r1_contig == r2_contig:
                    identity = (int(r1_match_bases) + int(r2_match_bases)) / (int(r1_alignment_block_length) + int(r2_alignment_block_length))
                    aligned_file.write("%s\t%s\t%s\t%s\t%s\n" %(r1_ID, r1_contig, r1_start_loc, r2_start_loc, identity))
                    r1_ID = "0"
                else:
                    r1_ID = "0"
    aligned_file.close()
    filtered_file = os.path.join(output_prefix, output_prefix + "_filtered_align.tsv")
    profile_file = os.path.join(output_prefix, output_prefix + "_profile.tsv")
    subprocess.run("Rscript rules/SGB_profile.r %s %s %s %s %s %s" %(cal_file,filtered_file,profile_file,genome_contig_length_file,identity,reads_num),shell=True)
    os.remove(cal_file)

def main():
    args = parser.parse_args()
    profile(args.minimap2_index,args.rmhost_reads1,args.rmhost_reads2,args.genome_contig_length,args.threads,args.identity,args.reads_num,args.output_prefix)

if __name__ == "__main__":
    main()
