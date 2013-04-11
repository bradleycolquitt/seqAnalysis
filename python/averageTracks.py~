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

# Combine subtraction and greater than zero check
#     to access arrays/lists only once  
def subtract_w_floor(one, two, floor_val):
    out_array = np.empty(len(one))
    out_array[:] = floor_val
    # Convert array to list. Substantial speedup.
    one = list(one)
    two = list(two)
    to_floor = True
    if floor_val < 0: to_floor = False
    for ind in range(len(one)):
        result = one[ind] - two[ind]
        if to_floor:
            if result > floor_val: out_array[ind] = result
        else:
            out_array[ind] = result
    return out_array
    
def run(one_track, two_track, out, to_floor):
    out_track_name = "_".join([one_track._v_name, "sub", two_track._v_name])
    print out_track_name
    test = checkIfNodeExists(out, out_track_name)
    if test: return
    one_track_chrs = one_track._f_listNodes()
    two_track_chrs = two_track._f_listNodes()
    #assert len(one_track_chrs) == len(two_track_chrs), "one chrs are %s, two chrs are %s" \
    #    % (one_track_chrs, two_track_chrs)
    floor_val = 0
    if not to_floor: floor_val = -1
    out_track_chr = []
    for chr in one_track._f_iterNodes():
        chr_name = chr._v_name
        print chr_name
        one_track_chr = one_track._f_getChild(chr_name)
        two_track_chr = two_track._f_getChild(chr_name)
        if to_floor:
            out_track_chr = subtract_w_floor(one_track_chr, two_track_chr, floor_val)
        else:
            #print "no floor"
            out_track_chr = one_track_chr[:] - two_track_chr[:]
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
    if atrack == "all":
        for one_track in one.iterNodes("/"):
            run(one_track, two_track, out, args.floor)
    else:
        one_track = one.getNode("/" + atrack)
        run(one_track, two_track, out, args.floor)

    out.flush()    
    out.close()
    
if __name__ == '__main__':
    main(sys.argv)
