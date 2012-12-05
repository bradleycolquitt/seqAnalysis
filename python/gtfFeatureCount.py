#!/usr/bin/env python

import sys
import argparse
import pdb

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', dest='gtf')
    parser.add_argument('-f', dest='feature')
    parser.add_argument('-k', dest='key', type=int, action="store")
    args = parser.parse_args()
    
    gtf = open(args.gtf)
    out_name = "_".join([args.gtf, args.feature])
    out = open(out_name, 'w')
    feature = args.feature
    key = args.key
    curr_key = ""
    count_dict = {}
    
    for line in gtf:
        pdb.set_trace()
        line = line.split()
        id = line[key]
        if id not in count_dict:
            count_dict[id] = 0
        if line[2] == feature:
            count_dict[id] = count_dict[id] + 1
    
    for record in iter(count_dict):
        out.write("\t".join([record.key, record.value]) + "\n")
    
    
if __name__ == '__main__':
    main(sys.argv)
    