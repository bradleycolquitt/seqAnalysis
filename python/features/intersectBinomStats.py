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
from scipy.stats import binom
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector




#GENOME_SIZE = 2654911517
GENOME_SIZE = 2000000000

def padjust(pvalues, method = "BH"):
    stats = importr('stats')
    p_adjust = stats.p_adjust(FloatVector(pvalues), method = method)
    return p_adjust

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
        n = int(line[1])    #reference_length
        k = int(line[4])    #count
        p = float(line[3]) / GENOME_SIZE  #feature_span
        pvalues[line[0]] = binom.pmf(k, n, p)
        
    qvalues = padjust(pvalues.values())
    
    output_path = "/".join([dir_path, "stats"])
    if not os.path.exists(output_path): os.mkdir(output_path)        
    out_file = open("/".join([output_path, "binom"]), 'w')
    out = "feature\tpvalue\tqvalue\n"
    out_file.write(out)
    
    for pvalue, qvalue in itertools.izip(pvalues.iteritems(), qvalues):
        #pvalues[key] = qvalue
        out = "{0}\t{1}\t{2}\n".format(pvalue[0], str(pvalue[1]), str(qvalue))
        out_file.write(out)
        
    
    
if __name__ == '__main__':
    main(sys.argv)
