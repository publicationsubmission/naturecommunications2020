#!/bin/sh

# for file in *.fasta.gz
# do
# 	metaphlan2.py $file ${file%.fasta.gz}_profile.txt --bowtie2out ${file%.fasta.gz}_bowtie2.txt --samout ${file%.fasta.gz}.sam.bz2 --input_type multifasta --nproc 36
# done

# E-Minus
# Sample1
# path="fq/"
# sample="E_minus_M1_S4"
# f1=$path$sample"_R1_001.fastq"
# f2=$path$sample"_R2_001.fastq"

# Sample2
# path="fq/"
# sample="E_minus_M2_S5"
# f1=$path$sample"_R1_001.fastq"
# f2=$path$sample"_R2_001.fastq"

# Sample3
path="fq/"
sample="E_minus_M3_S6"
f1=$path$sample"_R1_001.fastq"
f2=$path$sample"_R2_001.fastq"

# E-Plus
# Sample1
# path="fq/"
# sample="E_plus_M1_S1"
# f1=$path$sample"_R1_001.fastq"
# f2=$path$sample"_R2_001.fastq"

# Sample2
# path="fq/"
# sample="E_plus_M2_S2"
# f1=$path$sample"_R1_001.fastq"
# f2=$path$sample"_R2_001.fastq"

# Sample3
# path="fq/"
# sample="E_plus_M3_S3"
# f1=$path$sample"_R1_001.fastq"
# f2=$path$sample"_R2_001.fastq"

metaphlan2.py $f1,$f2 /tmp/$sample\_profile.txt --bowtie2out /tmp/$sample\_bowtie2.txt --samout /tmp/$sample.sam.bz2 --input_type fastq --nproc 36
