#!/usr/bin/env python

# Find bidirectional promoters
# Defined by TSSs within 5 kb of each other
# Specifically look for promoters with opposite strand genes

# Loop through gene BED file coordinate sorted
# Store previous gene TSS
# If current gene on opposite strand calculate TSS distance
# If less than 5kb, write both genes

import sys
import pdb

class bed_record:
    def __init__(self, string, gene):
        string_split = string.split()
        self.chrom = string_split[0]
        self.start = 0
        self.end = 0
        self.name = string_split[3]
        self.score = int(string_split[4])
        self.strand = string_split[5]
        
        if gene:
            if self.strand == "+":
                self.start = int(string_split[1])
                self.end = int(string_split[2])
            else:
                self.start = int(string_split[2])
                self.end = int(string_split[1])
        else:
            self.start = int(string_split[1])
            self.end = int(string_split[2])
            
        self.gene = gene
        
    ##def __init__(self):
    #    chrom = ""
    #    start = 0
    #    end = 0
    #    name = ""
    #    score = 0
    #    strand = ""
    
    def format(self):
        if self.gene:
            if self.strand == "-":
                start = self.end
                end = self.start
            else:
                start = self.start
                end = self.end
        return("\t".join([self.chrom, str(start), str(end), self.name, str(self.score), self.strand]) + "\n")

def main(argv):
    bed = open(argv[1])   
    bed_out = open(argv[1].split(".bed")[0] + "_bidir.bed", 'w')
    prev_record = bed_record("a\t0\t0\ta\t0\ta", True)
    
    for line in bed:
        curr_record = bed_record(line, True)
        if curr_record.chrom != prev_record.chrom:
            prev_record = curr_record
            continue
        else:
            if curr_record.strand != prev_record.strand:
                #pdb.set_trace()
                if (curr_record.start - prev_record.start) <= 5000:
                    bed_out.write(prev_record.format())
                    bed_out.write(curr_record.format())
            prev_record = curr_record    
                
                
        
if __name__ == '__main__':
    main(sys.argv)