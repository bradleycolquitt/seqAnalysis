#! /usr/bin/env python

import sys
import os
import re
import shutil
#import fileinput
#import glob
import string
import time
from operator import itemgetter
from multiprocessing import Process, Pool, Queue

global output_files
output_files = {}

global indices
indices = {}

def find_key(dic, val):
    """return the key of dictionary dic given the value"""    
    return [k for k, v in dic.iteritems() if re.search(v, val)]

## {{{ http://code.activestate.com/recipes/303060/ (r1)
def group(lst, n):
    """
    group([0,3,4,10,2,3], 2) => [(0,3), (4,10), (2,3)]
    
    Group a list into consecutive n-tuples. Incomplete tuples are
    discarded e.g.
    
    >>> group(range(10), 3)
    [(0, 1, 2), (3, 4, 5), (6, 7, 8)]
    """
    return zip(*[lst[i::n] for i in range(n)]) 
## end of http://code.activestate.com/recipes/303060/ }}}

def index_split(args):
    
    print args
    files = args[0]
    path = args[1]
    one_name = path + "/" + files[0]
    two_name = path + "/" + files[1]
    three_name = path + "/" + files[2]
    file_base = os.path.basename(one_name).split("_")[3]
    one = open(one_name)
    two = open(two_name)
    three = open(three_name)
    l1 = "a"
    while l1 != "":
        l1 = one.readline()
        l2 = two.readline().split()
        l3 = three.readline()
        if len(l2) > 10:
            # Check illumina PF
            if l2[10] == "1":
                # Check index read quality
                quality = list(l2[9])
                quality_convert = []
                for letter in quality:
                    quality_convert.append(ord(letter) - 65)
                below_thresh = 0
                for val in quality_convert:
                    if val < 20: below_thresh = below_thresh + 1
                if below_thresh >= 6:
                    continue
                """
                read1_quality = list(l1[9])
                read2_quality = list(l2[9])
                read1_quality_convert = []
                read2_quality_convert = []
                below_thresh1 = 0
                for letter in read1_quality:
                   if val < min_thresh: below_thresh = below_thresh + 1 
                """
                # Match index read to key
                index = find_key(indices, l2[8])
                if len(index) > 0: 
                    index = index[0]
                    index_path = "/".join([path, index])
                    if not os.path.exists(index_path): 
                        os.mkdir(index_path)
                    # Populate read file dict
                    if index not in output_files:
                        #print index
                        #print one_name, two_name, three_name
                        output_files[index] = [open("/".join([path, index, os.path.basename(one_name)]), 'a'),
                                               open("/".join([path, index, os.path.basename(two_name)]), 'a'),
                                               open("/".join([path, index, os.path.basename(three_name)]), 'a'),]                                                                                            
                    #Write to respective index read files
                    output_files[index][0].write(l1)
                    output_files[index][1].write("\t".join(l2) + "\n")
                    output_files[index][2].write(l3)
    for k, v in output_files.items():
        for i in v:
            i.close()
            i = ""
        v = ""    
    output_files.clear()

def main(argv):

    files = os.listdir(argv)
    test = re.compile("qseq")
    files = filter(test.search, files)

    ## Sort file names
    files_split = []
    for i in files:
        files_split.append(tuple(i.split("_")))
    files_split = sorted(files_split, key=itemgetter(3,2))

    files_join = []
    for i in files_split:
        files_join.append("_".join(i))
    files_group = group(files_join, 3)


    index = open("/seq/lib/illumina_index_sequences")
    for line in index:
        line = line.strip().split()
        indices[line[0]] = line[1]

    ## Split by indices
    pool = Pool(processes=10)
    start_time = time.time()
    #print files_group
    args = [(files, argv) for files in files_group]
    for arg in args:
        pool.apply_async(index_split, (arg,))
    #    index_split(files, argv)
    pool.close()
    pool.join()

    print time.time() - start_time
    
if __name__ == "__main__":
    main(sys.argv)



