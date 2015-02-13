#!/usr/bin/env python

import sys, re, pdb

def main(argv):
    
    infile = open(argv[1])
    outfile = open(argv[1] + ".txt", 'w')
    
    
    
    id = ""
    pen = ""
    left = ""
    right = ""
    flag_bad = False
    pen_re = re.compile("PRIMER_PAIR_0_PENALTY")
    
    for line in infile:
        if re.search("SEQUENCE_ID", line):
            id = line.split("=")[1].strip()
            continue
        
        if re.search("PRIMER_PAIR_NUM_RETURNED", line):
            num_ret = line.split("=")[1].strip()
            if float(num_ret) == 0:
                flag_bad = True
            else:
                flag_bad = False
            
        if pen_re.search(line):
            #pdb.set_trace()
            pen = line.split("=")[1].strip()
            continue
            #if float(pen) > 1: continue
            
        if re.search("PRIMER_LEFT_0_SEQUENCE", line):
            left = line.split("=")[1].strip()
        elif re.search("PRIMER_RIGHT_0_SEQUENCE", line):
            right = line.split("=")[1].strip()
        
        if line == "=\n" and not flag_bad:
            outfile.write("\t".join([id, left, right, pen]) + "\n")
            #left = ""
            #right = ""
        
if __name__ == '__main__':
    main(sys.argv)