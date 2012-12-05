#! /usr/bin/env python

import sys
import os
import re
import argparse
import shutil
import operator
import Bio.SeqIO as seq
import string
import time
import itertools
import pdb
from operator import itemgetter
from multiprocessing import Process, Pool, Queue


class indexer:
    def __init__(self, read1, read2, readInd, outpath, mismatches):
       # pdb.set_trace()
        self.read1 = os.path.abspath(read1)
        
        self.read2 = os.path.abspath(read2)
        self.readInd = os.path.abspath(readInd)
        self.mismatches = mismatches
        self.path_prefix = os.path.dirname(self.read1)
        self.read1f = open(read1)
        self.read2f = open(read2)
        self.readIndf = open(readInd)
        self.indices = {}
        self.outfiles = {}
        
        indicesf = open("/seq/lib/illumina_index_sequences")
        for index in indicesf:
            index = index.split()
            self.indices[int(index[0])] = index[1]
        
        #self.outpath = os.path.abspath(outpath)
        self.outpath = self.path_prefix + "/" + os.path.basename(read1).split("_")[0] + "_split"
        if not os.path.exists(self.outpath): os.makedirs(self.outpath)
        
        for index in self.indices.iterkeys():
            dir = "/".join([self.outpath, str(index)])
            if not os.path.exists(dir): os.makedirs(dir)
            self.outfiles[index] = [open("/".join([self.outpath, str(index), os.path.basename(read1)]), 'w'),
                                    open("/".join([self.outpath, str(index), os.path.basename(read2)]), 'w'),
                                    open("/".join([self.outpath, str(index), os.path.basename(readInd)]), 'w')]
        dir = "/".join([self.outpath, str(0)])
        if not os.path.exists(dir): os.makedirs(dir)    
        self.outfiles[0] = [open("/".join([self.outpath, str(0), os.path.basename(read1)]), 'w'),
                                    open("/".join([self.outpath, str(0), os.path.basename(read2)]), 'w'),
                                    open("/".join([self.outpath, str(0), os.path.basename(readInd)]), 'w')]   
        
            
    def split(self):
        rInd_seq = ""
        seq_dist1 = (0,)*12
        min1 = (0,)*2
        index = 0
        read_count = 0 
        for r1, r2, rInd in itertools.izip(seq.parse(self.read1f, "fastq"), seq.parse(self.read2f, "fastq"), seq.parse(self.readIndf, "fastq")):
            read_count = read_count + 1
            if read_count % 10000000 == 0: print>>sys.stderr, "Reads processed: " + str(read_count)
        
            rInd_seq = str(rInd.seq)        
            seq_dist1 = self.hammCompare(rInd_seq[:6])            
            min1 = min(enumerate(seq_dist1), key=itemgetter(1))
         
            index = 0
            if min1[1] <= self.mismatches: index = min1[0] + 1
            seq.write(r1, self.outfiles[index][0], "fastq")
            seq.write(r2, self.outfiles[index][1], "fastq")     
            seq.write(rInd, self.outfiles[index][2], "fastq")    

    def hammCompare(self, test_seq):
        return tuple(itertools.imap(hammD, self.indices.itervalues(), itertools.cycle([test_seq])))
 
#http://code.activestate.com/recipes/499304-hamming-distance/. This is the 'dalke' version.   
def hammD(str1, str2):
    #print str1, str2
    assert len(str1) == len(str2)
    ne = operator.ne
    return sum(itertools.imap(ne, str1, str2))   
    
def main(argv):
    parser = argparse.ArgumentParser(description="Split paired-end fastq files by index file.")
    parser.add_argument('-1', dest="read1")
    parser.add_argument('-2', dest="read2", required=False, default="none")
    parser.add_argument('-i', dest="readInd")
    parser.add_argument('-o', dest="outpath", default="none")
    parser.add_argument('-m', dest="mismatches", required=False, help="Number of mismatches tolerated in index sequence", default=1)
    args = parser.parse_args()
    #if not os.path.exists(args.outpath): os.mkdir(args.outpath)
    #error_log = open(args.outpath + "/error_log", 'a')
    #pdb.set_trace()
    i = indexer(args.read1, args.read2, args.readInd, args.outpath, args.mismatches)
    try:
        i.split()
    except:
        error_log.write(str(sys.exc_info()[0]) + "\n")
if __name__ == '__main__':
    main(sys.argv)


