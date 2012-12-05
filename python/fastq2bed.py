#!/usr/bin/env python

import os, sys
import argparse
import pysam
import bam2bed
import sam
from subprocess import Popen

def main(argv):
    sam_dir = "/media/storage2/data/sam/"
    bam_dir = "/media/storage2/data/bam/"
    bed_dir = "/media/storage2/data/bed/" 

    ## take Illmina fastq file
    parser = argparse.ArgumentParser(description="Map fastq files.")
    parser.add_argument('-i', '--input_fastq', nargs=2, metavar='fastq', required=True, 
                        dest='fastq_files', help='read1 and read3 of library')
    parser.add_argument('-n', '--index', dest='index', required=True, help='index number of library')
    parser.add_argument('-r', '--mate-inner-dist', metavar='mean', required=True, 
                        dest='mean', help='mean size of library (minus adaptors)')
    parser.add_argument('-s', '--mate-std-dev', metavar='sd', required=True, dest='sd', help='standard deviation of library')
    parser.add_argument('-g', '--GTF', metavar='gtf', dest='gtf', default='blank',
                        help='reference gtf file to map reads to')
    parser.add_argument('--stranded', action='store_true', required=False, dest='strand', default=False, help="indicate if library contains strand information")
    args = parser.parse_args()
    input1 = argv[1]
    input2 = argv[2]
    index = argv[3]
    input_prefix = input1.split("_")[0] + "_" + index
    print input_prefix
    #input_filtered = input_prefix + "_filt.fastq"
    #cmd_args = ['fastq_quality_filter', '-q', '20', '-p', '80', 
    #            '-i', input_fastq, '-o', input_filtered]
    #filt = Popen(cmd_args)
    #filt.wait()
    
    ## Run bowtie
    samfile = sam_dir + input_prefix + ".sam"
    if not os.path.exists(samfile):
        cmd_args = ['bowtie', '-S', '-m', '1', '-p', '6', '-I', '100', '-X', '1500', '--chunkmbs', '256', 
                    'mm9', '-1', input1, '-2', input2, samfile]
        print "Mapping with bowtie: " + " ".join(cmd_args[1:])
        bowtie = Popen(cmd_args)
        bowtie.wait()
    
    ## SAM to BAM
    bamfile = bam_dir + input_prefix + ".bam"
    if not os.path.exists(bamfile):
        sam.sam2bam(samfile, bamfile)
    ## Remove unmapped, duplicates, sort, index
    sam.proc(bam_dir, input_prefix)
    
    ## Output BED
    print "Writing BED..."
    bam = bam_dir + input_prefix + "_sort.bam"
    bed = bed_dir + input_prefix + ".bed"
    sam.bam2bed(bam, bed, pe = True)
    
if __name__ == "__main__":
    main(sys.argv)
