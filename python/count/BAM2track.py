#!/usr/bin/env python

import sys
import argparse
import pysam
import track_util
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

def run(bam, track_file_name, window_size):
    track = tb.openFile(track_file_name, "w")
    track_name = bam.split(".bam")[0]
    print(track_name)
    test = checkIfNodeExists(track, track_name)
    if test: return
    bam_file = pysam.Samfile(bam)
    chrs = bam_file.references
    vals = []
    chr_length = 0
    start = 0
    end = window_size
    i = 0
    for chr in chrs:
        chr_length = bam_file.lengths[bam_file.references.index(chr)]
        vals = np.zeros(chr_length / window_size)
        print chr
        #track_chr = track._f_getChild(chr_name)
        print ">> Counting..."
        while end < chr_length:
            vals[i] = float(bam_file.count(reference=chr, start=start, end=end)) / bam_file.mapped
            i += 1
            start = start + window_size
            end = end + window_size
        print ">> Creating array..."    
        track.createArray("/" + track_name, chr, vals)
 #       out_track_chr = out.getNode("/" + out_track_name, chr_name)
        track_util.setTrackAttributes(track, track_name + "/" + chr, 0, chr_length / window_size, chr, window_size)
#        def setTrackAttributes(file, node, start, stop, name, resolution):
#        for name in track_chr._v_attrs._f_list():
#            out.etNodeAttr("/" + "/".join([out_track_name, chr_name]), name, track_chr._v_attrs[name])
    print "Flushing..."
    track.flush()
    track.close()
    
def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', dest='bam', help='BAM file')
    parser.add_argument('-t', dest='track', help='track file')
    parser.add_argument('-w', dest='window_size', type=int, help='window_size')
    args = parser.parse_args()

    run(args.bam, args.track, args.window_size)
    
if __name__ == '__main__':
    main(sys.argv)
