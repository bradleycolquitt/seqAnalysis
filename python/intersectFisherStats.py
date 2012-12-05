#!/usr/bin/env python

# Read in intersection feature summary produced by intersectFeatureSet.py
# Compute fraction of feature set relative to whole genome
# Use this fraction as binomial probability
# K = number of peaks that intersect given set
# N = total number of peaks

import sys, os
import argparse
import itertools
import pdb
import rpy2.robjects as robjects
from scipy.stats import binom
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector




#GENOME_SIZE = 2654911517
GENOME_SIZE = 2000000000

stats = importr('stats')

def padjust(pvalues, method = "BH"):
    
    p_adjust = stats.p_adjust(FloatVector(pvalues), method = method)
    return p_adjust

def Rfisher(k, n, e):
    #print k,n,e
    fish = robjects.r['fisher.test']
    data = robjects.IntVector([k, n-k, e, n-e])
    fish_obj = fish(robjects.r['matrix'](data, ncol=2))
    #print fish_obj.names
    pvalue = fish_obj[0]
    return(pvalue)

def main(argv):
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='intersections')
    args = parser.parse_args()

    intersect_file = open(args.intersections)
    intersect_file.next()
    base_name = os.path.basename(args.intersections)
    dir_path = os.path.dirname(args.intersections)
    
    #pdb.set_trace()
    pvalues = {}
    for line in intersect_file:
        line = line.split()
        #print line[0]
        n = int(line[1])    #reference_length
        k = int(line[4])    #count
        p = float(line[3]) / GENOME_SIZE  #total size of feature space normalized by genome size
        #print p
        e = round(n * p) 
        #pvalues[line[0]] = binom.pmf(k, n, p)
        pvalues[line[0]] = Rfisher(k, n, e)
    print pvalues
    qvalues = padjust(pvalues.values())
    
    output_path = "/".join([dir_path, "stats"])
    if not os.path.exists(output_path): os.mkdir(output_path)        
    out_file = open("/".join([output_path, "binom"]), 'w')
    out = "feature\tpvalue\tqvalue\n"
    out_file.write(out)
    
    for pvalue, qvalue in itertools.izip(pvalues.iteritems(), qvalues):
        out = "{0}\t{1}\t{2}\n".format(pvalue[0], str(pvalue[1]), str(qvalue))
        out_file.write(out)
        
    
    
if __name__ == '__main__':
    main(sys.argv)
