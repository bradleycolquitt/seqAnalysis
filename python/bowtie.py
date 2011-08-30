#!/usr/bin/env python

import os, sys
import argparse
import pysam
import bam2bed
import sam
from subprocess import Popen

fastq_dir = "/media/storage2/data/fastq"
sam_dir = "/media/storage2/data/sam"
bam_dir = "/media/storage2/data/bam"
bed_dir = "/media/storage2/data/bed" 

class bowtie_class:
    def __init__(self, date, sample, single_end, index):
        self.date = date
        self.sample = sample.split("_")[1]
        self.single_end = single_end
        self.index = index
        self.input_prefix = index[1]
        self.input1 = ""
        self.input2 = ""
        if self.single_end:
            self.input1 = "/".join([fastq_dir, date, "".join([sample, '.fastq'])])
                                   
        else:
            self.input1 = "/".join([fastq_dir, date, sample, index[0], 
                                    "_".join([self.sample, '1.fastq'])])
            self.input2 = "/".join([fastq_dir, date, sample, index[0], 
                                    "_".join([self.sample, '3.fastq'])]) 
        self.samfile = "/".join([sam_dir, self.date, "".join([self.input_prefix, ".sam"])])
        #if os.path.exists(self.samfile):
        #    dec = raw_input("SAM file exists. Overwrite? [y/n]")
        #    if dec == "n": return
        #    elif dec == "y": os.remove(self.samfile)
        self.bamfile = "/".join([bam_dir, self.date, "".join([self.input_prefix, ".bam"])])
        sam_dir_date = "/".join([sam_dir, self.date])
        if not os.path.exists(sam_dir_date): os.mkdir(sam_dir_date)
        
        bam_dir_date = "/".join([bam_dir, self.date])
        if not os.path.exists(bam_dir_date): os.mkdir(bam_dir_date)

    def map(self):        
        if not os.path.exists(self.samfile):
            if not self.single_end:
                cmd_args = ['bowtie', '--phred33-quals', '-S', '-m', '1', '-p', '10', 
                            '-I', '100', '-X', '1500', '--chunkmbs', '256', 
                            'mm9', '-1', self.input1, '-2', self.input2, self.samfile]
            else:
                cmd_args = ['bowtie', '-S', '-m', '1', '-p', '10', 
                            '--chunkmbs', '256', 'mm9', self.input1, self.samfile]
            print "Mapping with bowtie: " + " ".join(cmd_args[1:])
            bowtie = Popen(cmd_args)
            bowtie.wait()
        else: print "Samfile exists. Skipping..."   
    def sam2bam(self):
        if not os.path.exists(self.bamfile):
            sam.sam2bam(self.samfile, self.bamfile)
        sam.proc(self.bamfile)
        
def bowtie(date, sample, single_end, index):
    bowtie_obj = bowtie_class(date, sample, single_end, index)
    bowtie_obj.map()
    bowtie_obj.sam2bam()
    
def main(argv):
    parser = argparse.ArgumentParser(description="Map fastq files.")
    parser.add_argument('-d', required=True, dest='date', help='sample date')
    parser.add_argument('-s', dest='sample', required=True, help='sample name')
    parser.add_argument('--single-end', action='store_true', type='bool', dest='single_end', default=False)
    parser.add_argument('-n', '--index', dest='index', required=True, help='index number and rename of library')
    args = parser.parse_args()
    
    bowtie(args.date, args.sample, args.single_end, args.index)

if __name__ == "__main__":
    main(sys.argv)
