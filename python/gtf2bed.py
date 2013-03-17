#!/usr/bin/env python

import sys
import pdb

def main(argv):
    
    infile = open(argv[1])
    outfile = open(argv[1].split(".gtf")[0] + ".bed", 'w')
    
    for line in infile:
        line = line.split('\t')
        chrom = line[0]
        start = line[3]
        end = line[4]
        strand = line[6]
        #pdb.set_trace()
       
        anno = line[8].split("; ")
        gene_id = anno[0].split("gene_id")[1].strip().split("\"")[1]
        gene_name = anno[4].split("gene_name")[1].strip().split("\"")[1]
        
        out = "\t".join([chrom, start, end, gene_name, gene_id, strand]) + "\n"
        outfile.write(out)
        
if __name__ == '__main__':
    main(sys.argv)