#!/bin/bash

conda activate mirna
####  Cutadapt to trim adapters 

find DIR -name "*fastq.gz" > raw_fastq ### replace DIR with the path to the raw fastq.gz 
raw_list="raw_fastq"
while read p
do
	sample=`echo $p | awk -F '\/' '{print $4,$5}' OFS="-" | sed s/_.*//g`  ####Adjust to extract the sample name
	mkdir trimmed/${sample}
	cutadapt -m 15 -u 3 -a AAAAAAAAAA -j 10 --trimmed-only $p > trimmed/${sample}/${sample}_trimmed.fastq.gz
done < $raw_list

find trimmed -name "*fastq" > trimmed_list


#### Align to known miRNA in Human
trim_list='trimmed_list'
while read p
do
	sample=`basename $p | sed s/_trim.*//g`
	mkdir Mirna_align/${sample}
	bowtie -q -v 0 -k 10 -p 10 --norc --best --strata -S ref/mirna_human $p Mirna_align/${sample}/${sample}.sam --un Mirna_align/${sample}/unalligned_${sample}.fastq
	samtools sort -o Mirna_align/${sample}/${sample}.bam Mirna_align/${sample}/${sample}.sam
	samtools index Mirna_align/${sample}/${sample}.bam
	samtools idxstats Mirna_align/${sample}/${sample}.bam | cut -f 1,3 > Mirna_align/${sample}/counts.tsv
	rm Mirna_align/${sample}/${sample}.sam
done < $trim_list

find Mirna_align -name "unalligned_*.fastq" > unalligned_mir_list

#### Align to whole genome in Human

un_list='unalligned_mir.list'
while read p
do
	sample=`echo $p | awk -F"\/" '{print $2}'`
	mkdir Genome_align/${sample}
	bowtie -q -n 1 -k 10 -p 10 --norc --best --strata -S ref/genome $p Genome_align/${sample}/${sample}.sam --un Genome_align/${sample}/unalligned_${sample}.fastq
	samtools sort Genome_align/${sample}/${sample}.sam > Genome_align/${sample}/${sample}.bam
	samtools index Genome_align/${sample}/${sample}.bam
	tagBam -i Genome_align/${sample}/${sample}.bam -files ref/target_mirna.bed -tag CM -names > Genome_align/${sample}/${sample}_tagged.bam
	rm Genome_align/${sample}/${sample}.sam
done < $un_list

find Genome_align -name "unalligned_*.fastq" > unalligned_gen_list

#### Align to MTB genome

un_list='unalligned_gen.list'
bowtie-build ref/mtb.fa ref/mtb


while read p
do
	sample=`echo $p | awk -F"\/" '{print $2}'`
	mkdir MTB_align/${sample}
	bowtie -q -v 0 -k 10 -p 10 --norc --best --strata -S ref/mtb $p MTB_align/${sample}/${sample}.sam --un MTB_align/${sample}/unalligned_${sample}.fastq
	samtools sort -o MTB_align/${sample}/${sample}.bam MTB_align/${sample}/${sample}.sam
	samtools index MTB_align/${sample}/${sample}.bam
	rm MTB_align/${sample}/${sample}.sam
done < $un_list


### Count with custom putative smRNA gtf (not shared)
find MTB_align -name "*.bam" > mtb_bam.list

featureCounts -t CDS -g gene_id -O -s 1 -M -a ref/mtb_all.gtf -o counts/mtb_count.tsv mtb_bam.list




