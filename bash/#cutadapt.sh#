#! /bin/bash

read1=(${1//.fastq/ })
read2=(${2//.fastq/ })

read1_tmp=${read1[0]}_tmp.fastq
read2_tmp=${read2[0]}_tmp.fastq

read1_trim=${read1[0]}_trimmed.fastq.gz
read2_trim=${read2[0]}_trimmed.fastq.gz

# filter index adaptor
cutadapt -a index=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --minimum-length 20 -o $read1_tmp -p $read2_tmp $1 $2 > ${read1[0]}_iadapt_trim.txt

# filter universal adaptor
cutadapt -a univ=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT --minimum-length 20 -o $read1_trim -p $read2_trim $read1_tmp $read2_tmp > ${read1[0]}_uadapt_trim.txt

rm $read1_tmp $read2_tmp
