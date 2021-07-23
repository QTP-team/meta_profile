## 1 Introduction
Accuracy and fast metagenomics classification using unique alignment.

Bowtie2 with the options: '--sensitive --end-to-end -k 2' is used to align reads to genomes.

Then alignments will be filtered according to the following criteria:
* Paired-reads align on the same contig;<br>
* unique alignment;<br>
* Alignment identity >= 0.95, Identity = 100 * MatchedBases / (MatchedBases + Substituions + Insertions + Deletions)<br>
* Genomes which supported by at least N (default: 100) paired-end reads are retained for further analysis.<br>

Finally, the relative abundance of each genome in the sample is calculated using the following formula:
rel_ab = (genome_alignment_num/genome_size) / sum(genome_alignment_num/genome_size) * 100

## 2 Usage
### 2.1 Use conda to install related dependencies
```
conda env create -n profile -f env.yaml
conda activate profile
```

### 2.2 Demo
The simulated metagenomics reads were generated using three genomes<br>
```
gzip -d 0.data/*.gz
iss generate --cpus 8 --draft 0.data/*.fa --n_reads 1M -m Hiseq --output test1_data
```

### 2.3 Index
Create bowtie2 index of reference genome<br>
```
mkdir 1.index
cat 0.data/*.fa > 1.index/merge_sp3.fna
bowtie2-build --threads 4 1.index/merge_sp3.fna 1.index/merge_sp3_bowtie2_index > 1.index/index.log
```

Obtain the corresponding relationship of genome_ID, contig_ID, contig_length.
```
python rules/genome_len.py -i 0.data/*.fa -o 1.index/genome_contig_length.txt
```

### 2.4 Run
Raw reads should be trimmed and removed the host contamination in the real study, omitted here.
```
python rules/profile_bowtie2.py -x 1.index/merge_sp3_bowtie2_index -r1 test1_data_R1.fastq -r2 test1_data_R2.fastq -contig 1.index/genome_contig_length.txt -t 8 -i 0.95 -n 100 -o test1
```

### 2.5 Output files
```test1_filtered_align.tsv```  Filtered alignments.<br>
```test1_profile.tsv```  Relative abundance of each genome in the sample.<br>
