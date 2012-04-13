#! /usr/bin/env python

import os, sys, re
import indexSplitter
import indexSplitter2
import argparse
import qseq
import periodToN
from subprocess import Popen
from multiprocessing import Pool

qseq_dir = '/media/storage3/data/qseq'
fastq_dir = '/media/storage3/data/fastq'

class index_class:
    def __init__(self, date, sample, se):
        self.date = date
        self.sample = sample
        self.se = se
        self.qseq_dir = "/".join([qseq_dir, date])
        self.fastq_dir = "/".join([fastq_dir, date])
        self.input_dir = self.qseq_dir + "/" + sample
        self.output_dir = "/".join([self.fastq_dir, sample])
        if not os.path.exists(self.fastq_dir): os.mkdir(self.fastq_dir)
        if not os.path.exists(self.output_dir): os.mkdir(self.output_dir)
    
    def split(self):
        if os.path.exists(self.input_dir + "/10"):
            pass
            #dec = raw_input("Split files exist. Split again [y/n]?")
            #if dec == "y": indexSplitter.main(self.input_dir)
        else:
            if self.se:
                indexSplitter2.main(self.input_dir)
            else:
                indexSplitter.main(self.input_dir)

    def convert(self):
        #args = [("/".join([self.input_dir, str(index)]), "/".join([self.output_dir, str(index)])) for index in range(0, 13)]
        #for arg in args:
        #    print arg
            #pool = Pool(processes = 10)
            #pool.apply_async(convert_worker, (arg,))
        #pool.close()
        #pool.join()
        #    qseq.convertToFastq(arg)    
        for index in range(0, 13):
            cmd_args = ['qseq2fastq', 
                    '-i', self.input_dir + "/" + str(index), 
                    '-o', "/".join([self.fastq_dir, self.sample, str(index)]),
                    '--preserve', '--threads', '6', '--filterpf', '--filterminqval', '20']
            """
            cmd_args = ['perl', '/home/user/src/perl/qseq2fastq.pl',
                        "/".join([self.input_dir, str(index), str(read)]),
                        "/".join([self.fastq_dir, self.sample, str(index), "".join([str(read), ".fastq"])])]
            """       
            p = Popen(cmd_args)
            p.wait()
        
    def pToN(self):
        for index in range(0, 13):
            skip = False
            files = os.listdir("/".join([self.output_dir, str(index)]))
            for file in files:
                if re.search("fastqc", file): skip = True
            if skip: continue    
            for file in files:
                if not re.search("count", file):
                    input = "/".join([self.output_dir, str(index), file])
                    print input
                    periodToN.main(['',input])
                    os.rename(input + "_Nconvert", input)
                    cmd_args = ['/seq/FastQC/fastqc', '--threads', '6', input]
                    p = Popen(cmd_args)
                    p.wait()
    
        
def index(date, sample, se):
    index_obj = index_class(date, sample, se)
    print "Splitting samples by indices..."
    index_obj.split()
    print "Converting to fastq..."
    index_obj.convert()
    print ". to N ..."
    index_obj.pToN()
    
    
    
def convert_worker(args):
    qseq.convertToFastq(args)
    
def main(argv):
    parser = argparse.ArgumentParser(description="Separate libraries by indices and convert to fastq.")
    parser.add_argument('-i', '--qseq', required=False, 
                        dest='qseq', help='qseq file')
    parser.add_argument('-d', required=True, dest="date", help='sample date')
    parser.add_argument('--single-end', dest="se", default=False, action="store_true", )
    args = parser.parse_args()
    
    if args.date and not args.qseq:
        qseqs = os.listdir(qseq_dir + "/" + args.date)
        for qseq in qseqs:
            print "Processing: " + qseq
            index(args.date, qseq, args.se)
    else:
        index(args.date, args.qseq, args.se)

if __name__ == "__main__":
    main(sys.argv)
