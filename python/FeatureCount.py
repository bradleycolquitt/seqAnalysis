#!/usr/bin/env python

import sys
import argparse
import pdb

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='input')
#    parser.add_argument('-f', dest='feature')
#    parser.add_argument('-k', dest='key', type=int, action="store")
    args = parser.parse_args()
    
    input = open(args.input)
    out_name = "_".join([args.input, "count"])
    out = open(out_name, 'w')
    #feature = args.feature
    #key = args.key
    #curr_key = ""
    count_dict = {}
    
    for line in input:
        #pdb.set_trace()
        line = line.split()
        name = line[3]
        name = name.split("_")
        id = name[1]
        if id not in count_dict:
            count_dict[id] = 0
        else:
            count_dict[id] = count_dict[id] + 1
    #pdb.set_trace()
    for record in count_dict.iteritems():
        out.write("\t".join(["_".join(["NM", record[0]]), str(record[1])]) + "\n")
    
    
if __name__ == '__main__':
    main(sys.argv)
    