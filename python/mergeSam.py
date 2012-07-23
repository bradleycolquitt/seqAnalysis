#!/usr/bin/env python

import sys
import argparse
from subprocess import Popen
def main(argv):
    
    parser = argparse.ArgumentParser(description="Map fastq files.")
    parser.add_argument('-i', required=True, nargs='+', dest='input', help='Sam files to be merged')
    parser.add_argument('-o', dest='output', required=True, help='Output name')
    parser.add_argument('--comment', required=False, help="Optional comment", default="")
    args = parser.parse_args()
    
    cmd_args = ['java', '-Xmx4g', '-jar', '/seq/picard/MergeSamFiles.jar',
                "=".join(["OUTPUT", args.output]),
                "ASSUME_SORTED=true", "SORT_ORDER=coordinate",
                "USE_THREADING=true",
                "=".join(["Comment", args.comment])] + \
                ["=".join(["INPUT", f]) for f in args.input]        
    p = Popen(cmd_args)
    p.wait()
    
if __name__ == '__main__':
    main(sys.argv)