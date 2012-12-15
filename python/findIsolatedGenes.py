#!/usr/bin/env python

# Find genes with TSSs distant from other genes
# Defined by TSSs at least 10 kb away from adjacent genes

# Loop through gene BED file coordinate sorted
# Store 2 previous genes
# Push gene back
# If middle gene is on plus strand, calculate distance from TSS to previous gene
# Else, calculate distance from TSS to next gene
# If distance is greater than 10 kb, write middle gene
# 

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

def calc_tss_distance(test, reference):
    val = 0
    if reference.strand == "+":
        val = reference.start - max(test.start, test.end)
    else:
        val = min(test.start, test.end) - reference.start
    return val

def main(argv):
    bed = open(argv[1])   
    bed_out = open(argv[1].split(".bed")[0] + "_dist10kb.bed", 'w')

    curr_record = ""
    prev_records = [bed_record("a\t0\t0\ta\t0\ta", True),
                    bed_record("chr1\t0\t0\ta\t0\ta", True)]
    tss_dist = 0
    
    initial_counter = 1
    for line in bed:
        curr_record = bed_record(line, True)
        
        if initial_counter > 0:
            prev_records[0] = prev_records[1]
            prev_records[1] = curr_record
            initial_counter -= 1
            continue
        
        # End of chromosome reached
        if curr_record.chrom != prev_records[1].chrom:
            if prev_records[1].strand == "+":
                tss_dist = calc_tss_distance(prev_records[0], prev_records[1])
            else:
                bed_out.write(prev_records[1].format())
            prev_records[0] = bed_record("a\t0\t0\ta\t0\ta", True)
            prev_records[1] = curr_record
        
        else:
            if prev_records[1].strand == "+":
                #pdb.set_trace()
                tss_dist = calc_tss_distance(prev_records[0], prev_records[1])
            else:
                tss_dist = calc_tss_distance(curr_record, prev_records[1])
                
            if (tss_dist) >= 10000:
                    bed_out.write(prev_records[1].format())
                
            prev_records[0] = prev_records[1]
            prev_records[1] = curr_record
                
                
        
if __name__ == '__main__':
    main(sys.argv)