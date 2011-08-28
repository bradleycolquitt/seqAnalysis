#!/usr/bin/env python

import sys, shutil, os, subprocess
from string import *

def main(argv):
    
    bedfile = open(argv[1], 'r')
    start_position = atoi(argv[2])
    end_position = atoi(argv[3])

    if os.path.exists("tmp"): 
        shutil.rmtree("tmp") 
    os.mkdir("tmp") 
    
    gene_file = open("tmp/genes", 'w')
    flank_file = open("tmp/flanks", 'w')

    for line in bedfile:
        sline = line.split()
        pos = atoi(sline[4])
        if pos >= start_position and pos <= end_position:
            gene_file.write(line)
        else:
            flank_file.write(line)
    gene_file.close()
    flank_file.close()
    
    #call intersectBed on gene_file, flank_file
    
    cmd_args = ['/usr/local/bin/intersectBed', '-a', 'tmp/flanks', '-b', 'tmp/genes', '-wa', '-v']
    nogene_file = open("tmp/nogene", 'w')
    p = subprocess.Popen(cmd_args,stdout=nogene_file)
    return_code = p.wait()
    nogene_file.close()

    gene_file = open("tmp/genes", 'r')
    nogene_file = open("tmp/nogene", 'r')
    final_out_name = argv[1] + "_nogene"
    final_out = open(final_out_name, 'w+')
    
    for line in gene_file:
        final_out.write(line)
    for line in nogene_file:
        final_out.write(line)

    gene_file.close()
    nogene_file.close()    
    final_out.close()
    #shutil.rmtree("tmp")
       
if __name__ == "__main__":
    main(sys.argv)
