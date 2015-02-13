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

def run(track, out, thresh):
    out_track_name = track._v_name + "_thresh" + str(thresh)
    print(out_track_name)
    test = checkIfNodeExists(out, out_track_name)
    if test: return
    track_chrs = track._f_listNodes()
    
    #assert len(track_chrs) == len(two_track_chrs), "one chrs are %s, two chrs are %s" \
    #    % (track_chrs, two_track_chrs)
    
    for chr in track._f_iterNodes():
       
        chr_name = chr._v_name
        print chr_name
        track_chr = track._f_getChild(chr_name)
        out_track_chr = track_chr[:]
        out_track_chr[out_track_chr[:] < thresh] = 0
        out.createArray("/" + out_track_name, chr_name, out_track_chr)
 #       out_track_chr = out.getNode("/" + out_track_name, chr_name)      
        for name in track_chr._v_attrs._f_list():
            out.setNodeAttr("/" + "/".join([out_track_name, chr_name]), name, track_chr._v_attrs[name])
    
def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', dest='track_file', help='track file')
    parser.add_argument('-n', dest='track_name', help='track name')
    parser.add_argument('-o', dest='out_name', help='out track file')
    parser.add_argument('--thresh', type=int, required=False, default=False)
    args = parser.parse_args()
    
    file = tb.openFile(args.track_file)

    out = tb.openFile(args.out_name, 'a')
    track = args.track_name
    
    
    if track == "all":
        for track in file.iterNodes("/"):
            run(track, out, args.thresh)
    else:
        track = one.getNode("/" + atrack)
        run(track, two_track, out, args.floor)

    out.flush()    
    out.close()
    
if __name__ == '__main__':
    main(sys.argv)
