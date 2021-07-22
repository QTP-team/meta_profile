# profile_bowtie2-v0.1
1 Introduction
It is a computational tool for profiling the composition of microbial communities (Bacteria, Archaea and Eukaryotes) from metagenomic shotgun sequencing data with species-level. The process includes: mapping of metagenomic sequencing data to a reference database, filtering of mapping results, and calculation of relative abundance of species.
  Advantage:
    Fast speed and high accuracy;
    Can customize the database;
    Not dependent on evolutionary relationships;
    Accurate estimation of species relative abundance.

1.1 Mapping
    Use bowtie2 for mapping, and filter the result that read unmapped, mate unmapped, and supplement alignment.

1.2 Filtering of mapping results
    Paired-reads match the same contig;
    unique alignment;
    Identity >= 0.95, Identity = 100 * MatchedBases / (MatchedBases + Substituions + Insertions + Deletions), MatchedBases = NM - Insertions - Deletions;
    Only keep SGBs with more than 100 reads on the comparison.

1.3 Calculation of relative abundance of species
    rel_ab = (SGB_reads_num/SGB_size) / sum(SGB_reads_num/SGB_size) * 100

2 Pre-requisites
2.1 Use conda to install related dependencies
    conda env create -n profile -f env.yaml
    conda activate profile
2.2 Demo data
The simulated metagenome were generated using three bacteria genomes
    sh 0.data/work.sh
Create bowtie2 index of reference genome
    sh 1.bowtie2_index/work.sh

3 Basic usage
    python profile.py -x bowtie2_index -r1 {rmhost_reads1} -r2 {rmhost_reads2} -contig {genome_contig_length_file} -ID {genome_ID_file} -t {threads} -i {identity} -n {reads_num} -o {output_prefix}

optional arguments:
  -h, --help            show this help message and exit
  -x BOWTIE2_INDEX, --bowtie2_index BOWTIE2_INDEX
                        The prefix of the Bowtie2 index of the reference genome
  -r1 RMHOST_READS1, --rmhost_reads1 RMHOST_READS1
                        Query input files are FASTQ .fq/.fastq
  -r2 RMHOST_READS2, --rmhost_reads2 RMHOST_READS2
                        Query input files are FASTQ .fq/.fastq
  -contig GENOME_CONTIG_LENGTH, --genome_contig_length GENOME_CONTIG_LENGTH
                        The file of the reference genomes corresponding to the length of contigs
  -ID GENOME_ID, --genome_ID GENOME_ID
                        The file containing the ID of the reference genome
  -t N, --threads N     Number of alignment threads to launch
  -i N, --identity N    Keep reads larger than N
  -n N, --reads_num N   Keep the genome larger than N reads
  -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                        The prefix name of the output result

3.1 Input files
rmhost_reads: Query input files are FASTQ .fq/.fastq.
bowtie2_index: The prefix of the Bowtie2 index of the reference genome.
genome_contig_length_file: The file of the reference genomes corresponding to the length of contigs.
genome_ID_file: The file containing the ID of the reference genome.

3.2 Run a single sample
sh work.sh

3.3 Output files
test1.bam: The result of the comparison of reads to the reference database.
test1_aligned.tsv: Keep paired-end reads aligned to the same contig and calculate the identity.
test1_profile.tsv: Filter the comparison results and calculate the relative abundance of species.
test1_summary.txt: The relative abundance of species in the statistical sample.