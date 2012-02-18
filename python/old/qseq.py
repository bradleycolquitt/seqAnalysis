#!/usr/bin/env python

from Bio import SeqIO

def convertToFastq(input, output):
    files = os.listdir(input)
    read1 = open("/".join([output, "1.fastq"]), 'w')
    read2 = open("/".join([output, "2.fastq"]), 'w')
    
    for file in files:
        indata = open(file)
        file_split = file.split("_")
        if int(file_split[2]) == 1:
            out = read1
        elif int(file_split[2]) == 3:
            out = read2
        for line in indata:
            line = line.strip().split()
            head = "@%s:%s:%s:%s:%s:%s#%s/%s PF=%s" % tuple(line[0:7] + [line[10]])
            read = line[8].replace(".", "N")
            qual = line[9]
            out.data("/n".join([head, read, head, qual]))
        indata.close()
    read1.close()
    read2.close()
    
#def illumina2sanger(input, output):

    