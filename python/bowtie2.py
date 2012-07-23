#!/usr/bin/env python

import os, sys
import argparse
import pysam
import bam2bed
import sam
import pdb
import datetime
from subprocess import Popen
from subprocess import PIPE

fastq_dir = "/media/storage3/data/fastq"
sam_dir = "/media/storage2/data/sam"
#sam_dir = "/home/user/data/sam"

bam_dir = "/media/storage2/data/bam"
bed_dir = "/media/storage2/data/bed" 

class bowtie_class:
    def __init__(self, date, sample, single_end, index, style):
        #pdb.set_trace()
        self.date = date
        self.sample = ""
        if style == "old":
            self.sample = sample.split("_")[1]
        elif style == "new":
            self.sample = sample
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
#        pdb.set_trace()
        #self.unmapped = "/".join([sam_date_dir, "".join([self.input_prefix, "_unaligned"])])
        self.bamfile = "/".join([bam_dir, self.date, "".join([self.input_prefix, ".bam"])])
        
        bam_date_dir = "/".join([bam_dir, self.date])
        if not os.path.exists(bam_date_dir): os.mkdir(bam_date_dir)
        fastq_date_log_dir = "/".join([fastq_date_dir, "logs"])
        #pdb.set_trace()
        if not os.path.exists(fastq_date_log_dir): os.mkdir(fastq_date_log_dir)
        #self.runlog = open("/".join([bam_date_log_dir, "".join([self.input_prefix, "_run_log"])]), 'w+')
        #pdb.set_trace()
        self.errorlog = open("/".join([fastq_date_log_dir, "".join([self.input_prefix, "_log"])]), 'a', 0)
        now = datetime.datetime.now()
        header = "[{0}/{1}/{2} {3}:{4}:{5}]\n".format(now.month, now.day,
                                                      now.year, now.hour,
                                                      now.minute, now.second)
        self.errorlog.write(header)
                          

    def map(self):
        #pdb.set_trace()
        run = True
        if os.path.exists(self.samfile):
            dec = raw_input("SAM exists. Overwrite? [y/n]")
            if dec == "n": run = False
        
        if run:    
            if not self.single_end:
                cmd_args = ['bowtie2',
                            '-p', '8',
                            '-I', '50', '-X', '1500',
                            #'--local', '--very-sensitive-local', '--mm',
                            '--end-to-end', '--mm',
                            '-x', 'mm9',
                            '-1', self.input1,
                            '-2', self.input2, 
                            '-S', self.samfile]
            else:
                cmd_args = ['bowtie', '-S', '-p', '6',
                            '--chunkmbs', '256', 'mm9', self.input1, self.samfile]
            self.errorlog.write(" ".join(cmd_args) + "\n")
            try:
                p1 = Popen(cmd_args, stderr=self.errorlog)
            #    print "Running bowtie: " + " ".join(cmd_args) 
                p1.wait()
            except:
                return
            #self.errorlog.close()
            
            #self.runlog.close()
    def sam2bam(self):
        #pdb.set_trace()
        if not os.path.exists(self.bamfile):
            try:
                sam.sam2bam(self.samfile, self.bamfile, self.errorlog)
            except IOError:
                self.errorlog.write("Can\'t open SAM for conversion.\n")
                return
            except:
                #pdb.set_trace()
                self.errorlog.write("SAM to BAM conversion failed\n")
                return
            else:
                self.errorlog.write("SAM to BAM completed successfully.\n")
                self.errorlog.write("Removing SAM...")
                os.remove(self.samfile)
            
            try:
                sam.proc([self.bamfile, "False", self.errorlog])
            except IOError:
                self.errorlog.write("BAM processing failed: IOError.\n")
            else:
                self.errorlog.write("BAM processing completed successfully.\n")
            
    def proc(self):
        ret = sam.proc([self.bamfile, "False"])

def bowtie(date, sample, single_end, index, style):
    #pdb.set_trace()
    bowtie_obj = bowtie_class(date, sample, single_end, index, style)
    bowtie_obj.map()
    #bowtie_obj.proc()
    bowtie_obj.sam2bam()
    
    
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
