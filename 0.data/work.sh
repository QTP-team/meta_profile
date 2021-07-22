gzip -d merge_sp3.fna.gz
iss generate --compress --cpus 8 -g merge_sp3.fna -b abundance_file.txt --n_reads 5M -m Hiseq --output test1_data
fastp -i test1_data_R1.fastq.gz -I test1_data_R2.fastq.gz -o test1_trim_R1.fastq.gz -O test1_trim_R2.fastq.gz -w 4 --length_required 70 -j test1.json -h test1.html 2> test1.log

samtools faidx merge_sp3.fna
grep '>' merge_sp3.fna | awk '{print $1}' | sed 's/>//g' > merge_sp3_ID.txt
while read -r line
do
    grep -w $line merge_sp3.fna.fai | awk '{print $1"\t"$2}' >> contig_length.txt
done < merge_sp3_ID.txt
paste merge_sp3_ID.txt contig_length.txt > genome_contig_length.txt
sed -i '1i SGB_ID\tContig_ID\tContig_length\t' genome_contig_length.txt
