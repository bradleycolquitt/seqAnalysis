#!/usr/bin/env python

import os
import sys
import argparse
import tables as tb
import numpy as np
import re
import pdb
import track_util as tutil

import genome.track
import genome.db


def create_carray(track, chrom, shape=None):
    atom = tb.UInt16Atom(dflt=0)

    zlib_filter = tb.Filters(complevel=1, complib="zlib")

    # create CArray for this chromosome
    if shape is None:
        shape = [chrom.length]
    else:
        shape = [shape]
    carray = track.h5f.createCArray(track.h5f.root, chrom.name,
                                    atom, shape, filters=zlib_filter)

    return carray

def run(in_track, out_track, gdb, scale, fnc):
    scale = int(scale)
    for chrom in gdb.get_chromosomes():
        out_length = chrom.length / scale
        out_carray = create_carray(out_track, chrom, out_length)
        inds = xrange(out_length)
        curr_array = in_track.get_nparray(chrom.name)
        max_val = out_length * scale
        out_array = fnc(curr_array[0:max_val].reshape(len(inds), scale), axis=1)
        out_carray[:] = out_array

    in_track.close()
    out_track.close()

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', dest='track', help='input track')
    parser.add_argument('-o', dest='out_name', help='out track file')
    parser.add_argument('-f', dest='function')
    parser.add_argument('-s', dest='scale')
    args = parser.parse_args()

    path = os.path.dirname(args.track)
    base = os.path.basename(args.track)


    out_fname = "/".join([path, args.function + args.scale, base])

    # open genome
    gdb = genome.db.GenomeDB()

    # open track
    in_track = gdb.open_track(args.track)

    if gdb.has_track(out_fname):
        gdb.delete_track(out_fname)
    out_track = gdb.create_track(out_fname)

    # shrink by specified factor and function
    if args.function == "mean":
        fnc = np.mean
    if args.function == "sum":
        fnc = np.sum
    run(in_track, out_track, gdb, args.scale, fnc)

if __name__ == '__main__':
    main(sys.argv)
