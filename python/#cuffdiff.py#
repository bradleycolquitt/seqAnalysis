#!/usr/bin/env python

import os, sys
import pysam
import bam2bed
import sam
import argparse
from subprocess import Popen

def main(argv):
    
    parser = argparse.ArgumentParser(description="Run cuffdiff on two tophat output files and a given gtf file")
    parser.add_argument('-a', '--first_file', metavar = 'bam1', required=True, dest='samp1', help='first bam file')
    parser.add_argument('-b', '--second_file', metavar = 'bam2', required=True, dest='samp2', help='second bam file')
    parser.add_argument('-g', '--gtf', dest='gtf', default='/seq/lib/refgene_name2.gtf', help='reference gtf file')
    parser.add_argument('-s', '--stranded', action='store_true', required=False, dest='strand', default=False, help="indicate if library contains strand information")
    args = parser.parse_args()
    
    ref_name = os.path.basename(args.gtf).split(".gtf")[0]
    samp1_name = os.path.basename(args.samp1).split(".bam")[0]
    samp2_name = os.path.basename(args.samp2).split(".bam")[0]
    library_type = 'fr-unstranded'
    if args.strand: library_type = 'fr-secondstrand'
    
    """
    ref_name = os.path.basename(argv[1]).split(".gtf")[0]
    samp1_name = os.path.basename(argv[2]).split(".bam")[0]
    samp2_name = os.path.basename(argv[3]).split(".bam")[0]
    """
    #cuffdiff_out =  "/media/storage2/data/rna/cuffdiff/" + samp1_name + "_" + samp2_name + "_" + ref_name
    #print cuffdiff_out
    ## Run cuffdiff
    """
    cmd_args = ['cuffdiff', '-o', cuffdiff_out, '-L', samp1_name + "," + samp2_name, 
                '-p', '8', '-N', '-v', argv[1], argv[2], argv[3]]
    """
    cmd_args = ['cuffdiff','-p', '8', '-L', samp1_name + ',' + samp2_name, '-N', 
                '--library-type', library_type, 
                args.gtf, args.samp1, args.samp2]
    print "Running cuffdiff: " + " ".join(cmd_args[1:])
    cuff = Popen(cmd_args)
    cuff.wait()

if __name__ == "__main__":
    main(sys.argv)
