#!/usr/bin/env python

import sys
import argparse
import re
import numpy as np

def changeStepSize(line, new_size):
    line = line.split()
    step_size = line[3].split("=")
    step_size[1] = new_size
    step_size = "=".join(step_size)
    line[3] = step_size
    return(" ".join(line) + "\n")
    
def main(argv):
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest="input", help="Input wig")
    parser.add_argument("-s", dest="new_step", help="New step size")
    args = parser.parse_args()

    in_wig = open(args.input)
    prefix = args.input.split(".wig")[0]
    out_name = prefix + "_step" + args.new_step + ".wig"
    out_wig = open(out_name, 'w')
    count = 0
    new_step_int = int(args.new_step)
    
    first_line = in_wig.readline()
    #print first_line
    first_line = first_line.split()
    step = first_line[3]
    step_size = step.split("=")
    agg_len = new_step_int / float(step_size[1])
    
    assert agg_len % 1 is not 0, "New step size isn't multiple of old"
    assert new_step_int > float(step_size[1]), "New step size isn't greater than old"
    
    step_size[1] = args.new_step
    step_size = "=".join(step_size)
    first_line[3] = step_size
    out_wig.write(" ".join(first_line) + "\n")
   
    agg = np.zeros(agg_len)
    
    head_flag = re.compile("f")
    
    for line in in_wig:
        if head_flag.match(line):
            #print line
            result = np.mean(agg[:count])
            out_wig.write(str(result) + "\n")
            count = 0
            out_wig.write(changeStepSize(line, args.new_step)) 
        else:
            if count < new_step_int:
                agg[count] = float(line.strip())
                count = count + 1
            else:
                result = np.mean(agg)
                if np.isnan(result): result = 0
                out_wig.write(str(result) + "\n")
                count = 0
if __name__ == '__main__':
    main(sys.argv)