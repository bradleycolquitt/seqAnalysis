#!/usr/bin/env python

import os, sys
import pysam
import bam2bed
import sam
from subprocess import Popen

def main(argv):
    #input_prefix = os.path.basename(os.getcwd())
    input_prefix = os.path.basename(argv[1]).split(".bam")[0]
     
    output_path = "/media/storage2/data/rna/cufflinks_strand/" + input_prefix
    
    if not os.path.exists(output_path): os.mkdir(output_path)
    
    ## Run cufflinks with given Tophat output and GTF file
    cmd_args = ['cufflinks', '-o', output_path, '-p', '10', '-N', 
                '--library-type', 'fr-secondstrand', argv[1]]
    print "Running cufflinks: " + " ".join(cmd_args[1:])

    cufflinks = Popen(cmd_args)
    cufflinks.wait()
  
if __name__ == "__main__":
    main(sys.argv)
