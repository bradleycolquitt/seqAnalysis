#!/usr/bin/env python

import os, sys
import pysam
import bam2bed
import sam
from subprocess import Popen

def main(argv):
    input_prefix = argv[1].split(".bam")[0]
    
    output_dir = "/media/Storage/user/data/rna/cufflinks" + 
    ## Run cufflinks with given Tophat output and GTF file
    cmd_args = ['cufflinks', '-S', '-m', '1', '-p', '8', 
                'mm9', argv[1], samfile]
    print "Mapping with bowtie: " + " ".join(cmd_args[1:])

    bowtie = Popen(cmd_args)
    bowtie.wait()
    ## Run cuffdiff
    
if __name__ == "__main__":
    main(sys.argv)
