#!/usr/bin/env python

import sys
import argparse
import re
import pdb
import numpy as np

def changeStepSize(line, new_size):
    line = line.split()
    step_size = line[3].split("=")
    step_size[1] = new_size
    step_size = "=".join(step_size)
    line[3] = step_size
    if len(line) == 5:
        span_size = line[4].split("=")
        span_size[1] = new_size
        span_size = "=".join(span_size)
        line[4] = span_size
    return(" ".join(line) + "\n")
    
def func_wrapper(func, *args):
    return(func(*args))
    
def main(argv):
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest="input", help="Input wig")
    parser.add_argument("-s", dest="new_step", help="New step size")
    parser.add_argument("-f", dest="fun", help="summary function")
    args = parser.parse_args()

    in_wig = open(args.input)
    prefix = args.input.split(".wig")[0]
    out_name = prefix + "_step" + args.new_step + "_" + args.fun + ".wig"
    out_wig = open(out_name, 'w')
    count = 0
    new_step_int = int(args.new_step)
    
    fun = ""
    if args.fun == "mean": fun = np.mean
    if args.fun == "var": fun = np.var
    
    first_line = in_wig.readline()
    #print first_line
    first_line = first_line.split()
    step = first_line[3]
    step_size = step.split("=")
    span_size = step_size
    if len(first_line) == 5:
        span = first_line[4]    
        span_size = span.split("=")

    agg_len = new_step_int / float(step_size[1])
    
    assert agg_len % 1 is not 0, "New step size isn't multiple of old"
    assert new_step_int > float(step_size[1]), "New step size isn't greater than old"
    
    step_size[1] = args.new_step
    step_size = "=".join(step_size)
    first_line[3] = step_size
    
    span_size[1] = args.new_step
    span_size = "=".join(span_size)
    if len(first_line) == 5:
        first_line[4] = span_size
    out_wig.write(" ".join(first_line) + "\n")
    
    agg = np.zeros(agg_len)
    
    head_flag = re.compile("f")
    
    for line in in_wig:
        #pdb.set_trace()
        if head_flag.match(line):
            #print line
            result = np.mean(agg[:count])
            out_wig.write(str(result) + "\n")
            count = 0
            out_wig.write(changeStepSize(line, args.new_step)) 
        else:
            if count < agg_len:
                agg[count] = float(line.strip())
                count = count + 1
            else:
                result = func_wrapper(fun, agg)
                if np.isnan(result): result = 0
                #if result > .1: pdb.set_trace()
                out_wig.write(str(result) + "\n")
                agg[0] = float(line.strip())
                count = 1
if __name__ == '__main__':
    main(sys.argv)