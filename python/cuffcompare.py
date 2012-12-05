#!/usr/bin/env python

import os, sys
import pysam
import bam2bed
import sam
from subprocess import Popen

def main(argv):
    input_prefix = os.path.basename(os.getcwd())
    #input_prefix = argv[1].split(".gtf")[0]
     
    output_path = "/media/storage2/data/rna/cuffcompare/" + input_prefix
    print output_path
    if not os.path.exists(output_path): os.mkdir(output_path)
    
    ## Run cufflinks with given Tophat output and GTF file
    cmd_args = ['cuffcompare', '-o', output_path, '-r', '/seq/lib/refgene_name2.gtf', '-R',
                argv[1]]
    print "Running cuffcompare: " + " ".join(cmd_args[1:])

    cuff = Popen(cmd_args)
    cuff.wait()
  
    ## Run cuffdiff
    #cuffdiff_out =  "/media/Storage/user/data/rna/cuffdiff/" + input_prefix
    #cmd_args = ['cuffdiff', '-o', cuffdiff_out, ]
if __name__ == "__main__":
    main(sys.argv)
