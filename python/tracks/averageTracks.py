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
    out_track_name = "_".join([one_track._v_name, "mean", two_track._v_name])
    print out_track_name
    test = checkIfNodeExists(out, out_track_name)
    if test: return
    one_track_chrs = one_track._f_listNodes()
    two_track_chrs = two_track._f_listNodes()
    
    for chr in one_track._f_iterNodes():
        chr_name = chr._v_name
        print chr_name
        one_track_chr = one_track._f_getChild(chr_name)
        two_track_chr = two_track._f_getChild(chr_name)
    
        # 2 x chr_length array
#        track_array = np.array([one_track_chr, two_track_chr])

        
        #out_track_chr = np.mean(track_array, axis=0)
        out_track_chr = (one_track_chr[:] + two_track_chr[:]) / 2
        out.createArray("/" + out_track_name, chr_name, out_track_chr)
        for name in one_track_chr._v_attrs._f_list():
            out.setNodeAttr("/" + "/".join([out_track_name, chr_name]), name, one_track_chr._v_attrs[name])
    
def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', dest='one_name', help='track file 1')
    parser.add_argument('-b', dest='two_name', help='track file 2')
    parser.add_argument('--atrack', dest='atrack', help="track name 1")
    parser.add_argument('--btrack', dest='btrack', help="track name 2")
    parser.add_argument('-o', dest='out_name', help='out track file')
    parser.add_argument('--floor', required=False, default=False)
    args = parser.parse_args()
    
    one = tb.openFile(args.one_name)
    two = tb.openFile(args.two_name)
    out = tb.openFile(args.out_name, 'a')
    atrack = args.atrack
    btrack = args.btrack
    two_track = two.getNode("/" + btrack)


    one_track = one.getNode("/" + atrack)
    run(one_track, two_track, out, args.floor)

    out.flush()    
    out.close()
    
if __name__ == '__main__':
    main(sys.argv)
