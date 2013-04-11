#!/usr/bin/env python

import sys
import argparse
import re
import pdb

def main(argv):
    ##
    parser = argparse.ArgumentParser()
    parser.add_argument('--bed')
    parser.add_argument('--tts')
    args = parser.parse_args()
    
    bed_file = open(args.bed)
    tts_file = open(args.tts)
    out_file = open(args.tts + ".bed", 'w')

    ## Create BED dict
    bed_dict = {}
    
    for line in bed_file:
        sline = line.split()
        bed_dict[sline[3]] = sline
        
    for line in tts_file:
        if not re.search("#", line):
            sline = line.split()
            tts_start = int(sline[1])
            tts_end = int(sline[2])
            try:
                bed_line = bed_dict[sline[0]]
            except KeyError:
                pdb.set_trace()
            
            out_line = "\t".join([bed_line[0],
                                  str(int(bed_line[1]) + tts_start),
                                  str(int(bed_line[1]) + tts_end),
                                  bed_line[3],
                                  "0",
                                  bed_line[5]]) + "\n"
                
            out_file.write(out_line)
        
    

if __name__ == '__main__':
    main(sys.argv)