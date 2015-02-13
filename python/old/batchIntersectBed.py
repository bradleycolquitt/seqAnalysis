#! /usr/bin/env python

import os, sys
import pysam
import bam2bed
import sam
from subprocess import Popen

def main(argv):
    a_bed = argv[1]
    b_bed = argv[2]
    
    a_bed_name = a_bed.split("/")[-1]
    a_bed_name = a_bed_name.split(".bed")[0]
    b_bed_name = b_bed_name = b_bed.split("/")[-1]
    out_bed_name = a_bed_name + "_inter_" + b_bed_name
    out_bed = open(out_bed_name, 'w')
    cmd_args = ['intersectBed', '-a', a_bed, '-b', b_bed, '-c', '-wa']
    
    inter = Popen(cmd_args, stdout = out_bed)
    print "Intersecting " + a_bed + " with " + b_bed
    inter.wait()
    
if __name__ == "__main__":
    main(sys.argv)
