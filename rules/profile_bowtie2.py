import argparse
import pysam
import re
import subprocess
import os

parser = argparse.ArgumentParser(
    prog = "profile_bowtie2.py",
    usage = "python rules/profile_bowtie2.py -x 1.index/merge_sp3_bowtie2_index -r1 test1_data_R1.fastq -r2 test1_data_R2.fastq -contig 1.index/genome_contig_length.txt -t 8 -i 0.95 -n 100 -o test1")

arg = parser.add_argument

arg('-x', '--bowtie2_index', type=str, help = 'The prefix of the Bowtie2 index of the reference genome')
arg('-r1', '--rmhost_reads1', type=str, help = 'Query input files are FASTQ .fq/.fastq')
arg('-r2', '--rmhost_reads2', type=str, help = 'Query input files are FASTQ .fq/.fastq')
arg('-contig', '--genome_contig_length', type=str, help = 'The file of the reference genomes corresponding to the length of contigs')
arg('-t', '--threads', metavar = 'N', type=int, default=4, help = 'Threads for bowtie2')
arg('-id', '--identity', metavar = 'N', type=float, default=0.95, help = 'Alignments which identity larger than N will be retained')
arg('-n', '--reads_num', metavar = 'N', type=int, default=100, help = 'Genomes which supported by at least N paired-end reads will be retained')
arg('-o', '--output_prefix', type=str, help = 'The prefix of the output result')


def profile(index,rmhost1,rmhost2,genome_contig_length_file,threads,identity,reads_num,output_prefix):
    os.mkdir(output_prefix)
    map_file = os.path.join(output_prefix, output_prefix + ".bam")
    bam_file = open(map_file, "w")
    cal_file = os.path.join(output_prefix, output_prefix + "_tmp.tsv")
    aligned_file = open(cal_file, "w")
    subprocess.run("bowtie2 -p %s --sensitive --end-to-end -k 2 -x %s -1 %s -2 %s | samtools view -F 2060 -b" %(threads,index,rmhost1,rmhost2),shell=True,encoding="utf-8",stdout=bam_file)
    bam_file.close()
    aligned_file.write("Read_ID\tContig_ID\tR1_start_loc\tR2_start_loc\tIdentity\n")
    with pysam.AlignmentFile(map_file, "rb") as f:
        n = 0
        for line in f:
            n += 1
            if (n % 2) != 0:
                r1_ID = line.query_name
                r1_start_loc = line.reference_start
                r1_CIGAR = line.cigarstring
                r1_contig = line.reference_name
                r1_NM = line.get_tag('NM')
                r1_match_mismatch_bases = sum([ int(x) for x in re.findall(r'(\d+)M', r1_CIGAR) ])
                r1_InDel_bases = sum([ int(x) for x in re.findall(r'(\d+)[ID]', r1_CIGAR) ])
                r1_match_bases = int(r1_match_mismatch_bases) - int(r1_NM) + int(r1_InDel_bases)
            else:
                r2_ID = line.query_name
                r2_start_loc = line.reference_start
                r2_CIGAR = line.cigarstring
                r2_contig = line.reference_name
                r2_NM = line.get_tag('NM')
                if r1_ID == r2_ID and r1_contig == r2_contig:
                    r2_match_mismatch_bases = sum([ int(x) for x in re.findall(r'(\d+)M', r2_CIGAR) ])
                    r2_InDel_bases = sum([ int(x) for x in re.findall(r'(\d+)[ID]', r2_CIGAR) ])
                    r2_match_bases = int(r2_match_mismatch_bases) - int(r2_NM) + int(r2_InDel_bases)
                    r1_alignment_block_length = sum([ int(x) for x in re.findall(r'(\d+)[MDI]', r1_CIGAR) ])
                    r2_alignment_block_length = sum([ int(x) for x in re.findall(r'(\d+)[MDI]', r2_CIGAR) ])
                    identity = (r1_match_bases + r2_match_bases) / (int(r1_alignment_block_length) + int(r2_alignment_block_length))
                    aligned_file.write("%s\t%s\t%s\t%s\t%s\n" %(r1_ID, r1_contig, r1_start_loc, r2_start_loc, identity))
    aligned_file.close()
    filtered_file = os.path.join(output_prefix, output_prefix + "_filtered_align.tsv")
    profile_file = os.path.join(output_prefix, output_prefix + "_profile.tsv")
    subprocess.run("Rscript rules/SGB_profile.r %s %s %s %s %s %s" %(cal_file,filtered_file,profile_file,genome_contig_length_file,identity,reads_num),shell=True)
    os.remove(cal_file)

def main():
    args = parser.parse_args()
    profile(args.bowtie2_index,args.rmhost_reads1,args.rmhost_reads2,args.genome_contig_length,args.threads,args.identity,args.reads_num,args.output_prefix)

if __name__ == "__main__":
    main()
