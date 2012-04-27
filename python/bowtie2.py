#!/usr/bin/env python

import os, sys
import argparse
import pysam
import bam2bed
import sam
import pdb
from subprocess import Popen

fastq_dir = "/media/storage3/data/fastq"
#sam_dir = "/media/storage2/data/sam"
sam_dir = "/home/user/data/sam"
bam_dir = "/media/storage2/data/bam"
bed_dir = "/media/storage2/data/bed" 

class bowtie_class:
    def __init__(self, date, sample, single_end, index):
        #pdb.set_trace()
        self.date = date
        self.sample = sample.split("_")[1]
        #self.sample = sample
        self.single_end = single_end
        self.index = index
        self.input_prefix = index[1]
        self.input1 = ""
        self.input2 = ""
        fastq_date_dir = "/".join([fastq_dir, date])
        
        if single_end:
            self.input1 = "/".join([fastq_dir, date, "".join([sample, '.fastq'])])                      
        else:
            self.input1 = "/".join([fastq_dir, date, sample, index[0], 
                                    "_".join([self.sample, '1.fastq'])])
            self.input2 = "/".join([fastq_dir, date, sample, index[0], 
                                    "_".join([self.sample, '3.fastq'])])
            print self.input1
            #self.input1 = "/".join([fastq_dir, date, sample, index[0], '1.fastq'])
            #self.input2 = "/".join([fastq_dir, date, sample, index[0], '2.fastq'])
        sam_date_dir = "/".join([sam_dir, self.date])
        if not os.path.exists(sam_date_dir): os.mkdir(sam_date_dir)
        self.samfile = "/".join([sam_date_dir, "".join([self.input_prefix, ".sam"])])
        sam_date_log_dir = "/".join([sam_date_dir, "log"])
        
        if not os.path.exists(sam_date_log_dir): os.mkdir(sam_date_log_dir)
        self.errorlog = "/".join([sam_date_log_dir, "".join([self.input_prefix, "_log"])])
        #self.unmapped = "/".join([sam_date_dir, "".join([self.input_prefix, "_unaligned"])])
        self.bamfile = "/".join([bam_dir, self.date, "".join([self.input_prefix, ".bam"])])
        
        bam_dir_date = "/".join([bam_dir, self.date])
        if not os.path.exists(bam_dir_date): os.mkdir(bam_dir_date)

    def map(self):        
        if not os.path.exists(self.samfile):
            if not self.single_end:
                cmd_args = ['bowtie2',
                            '-p', '6', 
                            '-I', '50', '-X', '1500',
                            '--end-to-end',
                            '-x', 'mm9',
                            '-1', self.input1,
                            '-2', self.input2,
                            '|', 'samtools', 'view', '-Sb', '-o', self.bamfile, '-']
                cmd_args2 = ['samtools', 'view', '-u', '-F', '4', self.bamfile, '|',
                             'samtools', 'sort', '-', self.input_prefix + "_sort"]
                            
            else:
                cmd_args = ['bowtie', '-S', '-p', '8',
                            '--chunkmbs', '256', 'mm9', self.input1, self.samfile]
            print "Mapping with bowtie: " + " ".join(cmd_args[1:])
            errorlog = open(self.errorlog, 'w')
            bowtie = Popen(cmd_args, stderr=errorlog)
            bowtie.wait()
            samtools = Popen(cmd_args, stderr=errorlog)
            #errorlog.close()
            
    def sam2bam(self):
        if not os.path.exists(self.bamfile):
            sam.sam2bam(self.samfile, self.bamfile)
            ret = sam.proc([self.bamfile, "False"])
            if ret == 0: os.remove(self.samfile)
            #print "sam2bam"
            #sam.proc_sam([self.samfile, "True"])
            #samfile = self.samfile.split(".sam")[0]
            #sam.sam2bam(samfile + "_sort.sam", self.bamfile)
        
def bowtie(date, sample, single_end, index):
    bowtie_obj = bowtie_class(date, sample, single_end, index)
    #pdb.set_trace()
    bowtie_obj.map()
    #bowtie_obj.sam2bam()
    
    
def main(argv):
    parser = argparse.ArgumentParser(description="Map fastq files.")
    parser.add_argument('-d', required=True, dest='date', help='sample date')
    parser.add_argument('-s', dest='sample', required=True, help='sample name')
    parser.add_argument('-1', dest='fastq1')
    parser.add_argument('-2', dest='fastq2')
    parser.add_argument('--single-end', action='store_true', dest='single_end', default=False)
    parser.add_argument('-n', '--index', dest='index', required=True, help='index number of library')
    args = parser.parse_args()
    index_split = args.index.split("-")
    bowtie(args.date, args.sample, args.single_end, index_split)

if __name__ == "__main__":
    main(sys.argv)
