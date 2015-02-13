#!/usr/bin/env python

import sys, os
import argparse
import math
import tables as tb
from multiprocessing import Pool

def worker(file):
    print file
    h5 = tb.openFile(file, 'r')
    copy_name = file.split(".trk")[0] + "_copy.trk"
    h5.copyFile(copy_name)
    h5.close()
    h5 = tb.openFile(copy_name, 'a')

    samples = h5.iterNodes("/")
    for sample in samples:
        chrs = sample._f_iterNodes()
        for chr in chrs:
            for ind in xrange(len(chr)):
                if math.isnan(chr[ind]):
                    chr[ind] = 0
    h5.close()   

def main(argv):    
    files = os.listdir(argv[1])
    print files
    #for file in files:
    #    worker(file)
    pool = Pool(processes=3)
    [pool.apply_async(worker, (file,)) for file in files]
    pool.close()
    pool.join()

if __name__ == '__main__':
    main(sys.argv)