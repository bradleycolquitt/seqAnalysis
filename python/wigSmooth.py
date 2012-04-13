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
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest="input", help="Input wig")
    parser.add_argument("-s", dest="window_size", help="Smoothing window size")
    args = parser.parse_args()

    in_wig = open(args.input)
    prefix = args.input.split(".wig")[0]
    out_name = prefix + "_smooth" + args.window_size + ".wig"
    out_wig = open(out_name, 'w')
    count = 0
    window_size = int(args.window_size)
    
    ##Skip first header. If header encountered, reset window
    ##Save number of values specified in window size
    ##Mean values and write
    
    
    first_line = in_wig.readline()
    out_wig.write(first_line)
    ##print first_line
    #first_line = first_line.split()
    #step = first_line[3]
    #step_size = step.split("=")
    #agg_len = new_step_int / float(step_size[1])
    #
    #assert agg_len % 1 is not 0, "New step size isn't multiple of old"
    #assert new_step_int > float(step_size[1]), "New step size isn't greater than old"
    #
    #step_size[1] = args.new_step
    #step_size = "=".join(step_size)
    #first_line[3] = step_size
    #out_wig.write(" ".join(first_line) + "\n")
    #
    #agg = np.zeros(agg_len)
    
    head_flag = re.compile("f")
    window_arr = np.zeros(window_size)
    initial_fill_count = 0
    smoothed_val = 0
    for line in in_wig:
        if head_flag.match(line):
            print line
            smoothed_val = np.mean(window_arr)
            out_wig.write(str(smoothed_val) + "\n")
            out_wig.write(line)
            initial_fill_count = 0
            #out_wig.write(changeStepSize(line, args.new_step)) 
        else:
            val = float(line.strip())
            if initial_fill_count < window_size:
                window_arr[initial_fill_count] = val
                initial_fill_count = initial_fill_count + 1
            else:
                smoothed_val = np.mean(window_arr)
                out_wig.write(str(smoothed_val) + "\n")
                window_arr = np.roll(window_arr, -1)
                window_arr[window_size - 1] = val
                
if __name__ == '__main__':
    main(sys.argv)