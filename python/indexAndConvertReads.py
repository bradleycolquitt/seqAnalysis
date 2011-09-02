#! /usr/bin/env python

import os, sys
import indexSplitter
import argparse
from subprocess import Popen

qseq_dir = '/media/storage2/data/qseq'
fastq_dir = '/media/storage2/data/fastq'

class index_class:
    def __init__(self, date, sample):
        self.date = date
        self.sample = sample
        self.qseq_dir = "/".join([qseq_dir, date])
        self.fastq_dir = "/".join([fastq_dir, date])
        self.input_dir = self.qseq_dir + "/" + sample
        self.output_dir = "/".join([self.fastq_dir, sample])

        if not os.path.exists(self.fastq_dir): os.mkdir(self.fastq_dir)
        if not os.path.exists(self.output_dir): os.mkdir(self.output_dir)
    
    def split(self):
        if os.path.exists(self.input_dir + "/1"):
            dec = raw_input("Split files exist. Split again [y/n]?")
            if dec == "y": indexSplitter.main(self.input_dir)
        else:
            indexSplitter.main(self.input_dir)

    def convert(self):
        for index in range(0,7):
            split_index_dir = self.input_dir + "/" + str(index)
            if os.path.exists(split_index_dir):
                index_dir = self.output_dir + "/" + str(index)
                if not os.path.exists(index_dir): os.mkdir(index_dir)
                else: continue
                cmd_args = ['qseq2fastq', 
                        '-i', self.input_dir + "/" + str(index), 
                        '-o', "/".join([self.fastq_dir, self.sample, str(index)]),
                        '--preserve', '--threads', '10']
                p = Popen(cmd_args)
                p.wait()
     
def index(date, sample):
    index_obj = index_class(date, sample)
    print "Splitting samples by indices..."
    index_obj.split()
    print "Converting to fastq..."
    index_obj.convert()

def main(argv):
    parser = argparse.ArgumentParser(description="Separate libraries by indices and convert to fastq.")
    parser.add_argument('-i', '--qseq', required=False, 
                        dest='qseq', help='qseq file')
    parser.add_argument('-d', required=True, dest="date", help='sample date')
    args = parser.parse_args()
    
    if args.date and not args.qseq:
        qseqs = os.listdir(qseq_dir + "/" + args.date)
        for qseq in qseqs:
            print "Processing: " + qseq
            index(args.date, qseq)
    else:
        index(args.date, args.qseq)

if __name__ == "__main__":
    main(sys.argv)
