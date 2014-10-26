#! /bin/bash

cutadapt_run () {
    declare -a read1=(`ls | grep fastq.gz | grep R1`)
    declare -a read2=(`ls | grep fastq.gz | grep R2`)

    prefix1=(${read1[$1]//.fastq.gz/ })
    prefix2=(${read2[$1]//.fastq.gz/ })    
 
    echo $prefix1
    echo $prefix2
    read1_tmp=trimmed/${prefix1}_tmp.fastq
    read2_tmp=trimmed/${prefix2}_tmp.fastq
 
    read1_trim=trimmed/${prefix1}_trimmed.fastq.gz
    read2_trim=trimmed/${prefix2}_trimmed.fastq.gz

    # filter index adaptor
    echo "Filter index adaptor"
    cutadapt -a index=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --minimum-length 20 -o $read1_tmp -p $read2_tmp ${read1[$1]} ${read2[$1]} > trimmed/${prefix1}_iadapt_trim.txt

    # filter universal adaptor
    echo "Filter universal adaptor"
    cutadapt -a univ=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT --minimum-length 20 -o $read1_trim -p $read2_trim $read1_tmp $read2_tmp > trimmed/${prefix1}_uadapt_trim.txt

    rm $read1_tmp $read2_tmp
}

declare -a read1=(`ls | grep fastq.gz | grep R1`)
array_length=$((${#read1[@]} - 1))
echo ${read1[*]}

# Argument 1: number of processes to use
nproc=$1

mkdir trimmed
export -f cutadapt_run
seq 0 $array_length | parallel --progress --verbose -I index -j $nproc "cutadapt_run index"
