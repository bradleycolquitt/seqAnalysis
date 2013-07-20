#!/usr/bin/env python

import sys
import argparse
import re
import numpy as np

def smooth(line, new_size):
    line = line.split()
    step_size = line[3].split("=")
    step_size[1] = new_size
    step_size = "=".join(step_size)
    line[3] = step_size
    return(" ".join(line) + "\n")
    
def main(argv):
    
    #parser = argparse.ArgumentParser()
    #parser.add_argument("-i", dest="input", help="Input wig")
    #parser.add_argument("-o", dest="output")
    #args = parser.parse_args()

    in_wig = open(argv[1])
    prefix = argv[1].split(".wig")[0]
    out_name = prefix + "_var.wig"
    out_wig = open(out_name, 'w')
    count = 0
    
    
    ##Skip first header. If header encountered, reset window
    ##Save number of values specified in window size
    
    first_line = in_wig.readline()
    first_line_split = first_line.split()
    window_size = int(first_line_split[3].split("=")[1])
    first_line = " ".join(["variableStep"] + first_line_split[1:]) + "\n"
    out_wig.write(first_line)

    head_flag = re.compile("f")
    window_arr = np.zeros(window_size)
    initial_fill_count = 0
    pos = 1

    for line in in_wig:            
        if head_flag.match(line):
            print line
            line_split = line.split()
            out_wig.write(" ".join(["variableStep"] + line_split[1:]))
            pos = 0
        else:
            val = float(line.strip())   
            out_wig.write("\t".join([str(pos), str(val)]) + "\n")
            pos = pos + window_size    
                
if __name__ == '__main__':
    main(sys.argv)
