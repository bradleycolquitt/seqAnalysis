#!/usr/bin/env python


## Originally written to remove first and last exons within gene group
## Infile is 
import sys
import pdb

def main(argv):
    infile = open(argv[1], 'r')
    outfile = open(argv[2], 'w')
    first = True
    trans_store = ""
    line_store = ""
    for line in infile:
        #pdb.set_trace()
        sline = line.split()
        #anno = sline[8].split("; ")
        curr_trans = sline[11]
        
        #if first:
        #    trans_store = curr_trans
            #line_store = line
        #    first = False
        #    continue
        
        if curr_trans != trans_store:    
            trans_store = curr_trans
            line_store = ""
            #first = True
            continue
        else:
            if line_store != "":
                outfile.write(line_store)
            line_store = line
        
        
if __name__ == '__main__':
    main(sys.argv)
