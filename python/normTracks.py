#!/usr/bin/env python

import sys
import argparse
import tables as tb
import numpy as np
import itertools
import re

def checkIfNodeExists(out, name):
    nodes = out.listNodes("/")
    try:
        out.getNode("/" + name)
    except tb.NoSuchNodeError:
        pass
    else:
        dec = raw_input("Node exists. Overwrite [y/n]? ")
        if dec == "n": return True
        elif dec == "y":
            out.removeNode("/" + name, recursive=True)
    out.createGroup("/", name)
    return False
    
def run(one_track, two_track, out, to_floor):
    out_trace_name = one_track._v_name + "_norm"
    print out_track_name
    test = checkIfNodeExists(out, out_track_name)
    if test: return
    
    one_track_chrs = one_track._f_listNodes()
   
    value_total = 0 
    out_track_chr = []
    for chr in one_track._f_iterNodes():
        chr_name = chr._v_name
        print chr_name
        
        one_track_chr = one_track._f_getChild(chr_name)
        value_total += sum(one_track_chr[:])
    
    for chr in track._f_iterNodes():
        chr_name = chr._v_name
        print chr_name
        one_track_chr = one_track._f_getChild(chr_name)
        
        out_track_chr = 1E6 * one_track_chr / value_total
        out.createArray("/" + out_track_name, chr_name, out_track_chr)
        for name in one_track_chr._v_attrs._f_list():
            out.setNodeAttr("/" + "/".join([out_track_name, chr_name]), name, one_track_chr._v_attrs[name])
    
def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', dest='one_name', help='track file 1')
    parser.add_argument('--atrack', dest='atrack', help="track name 1")
    parser.add_argument('-o', dest='out_name', help='out track file')
    args = parser.parse_args()
    
    one = tb.openFile(args.one_name)
    out = tb.openFile(args.out_name, 'a')
    atrack = args.atrack

    if atrack == "all":
        for one_track in one.iterNodes("/"):
            run(one_track, out)
    else:
        one_track = one.getNode("/" + atrack)
        run(one_track, out)

    out.flush()    
    out.close()
    
if __name__ == '__main__':
    main(sys.argv)
