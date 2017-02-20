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

def run(track_fname, gdb, norm_type):
    in_track = gdb.open_track(track_fname)
    path = os.path.dirname(track_fname)
    base = os.path.basename(track_fname)
    out_fname = "/".join([path, norm_type, base])
    gdb.delete_track(out_fname)
    out_track = gdb.create_track(out_fname)

    value_total = 0

    print "Counting reads..."
    for chrom in gdb.get_chromosomes():
        value_total += np.sum(in_track.get_nparray(chrom.name))

    norm_factor = 1
    if norm_type == 'rpm':
        norm_factor = float(value_total) / 1E6

    print "Norming reads..."
    for chrom in gdb.get_chromosomes():
        out_carray = create_carray(out_track, chrom)
        curr_array = in_track.get_nparray(chrom.name) / norm_factor
        out_carray[:] = curr_array

    in_track.close()
    out_track.close()

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', dest='track', help='input track')
    parser.add_argument('--assembly', help="assembly to use", default='lonStrDom1')
    parser.add_argument('-n', dest='norm_type', choices=['rpm'], default='rpm')
    args = parser.parse_args()

    # open genome
    gdb = genome.db.GenomeDB()
    run(args.track, gdb, args.norm_type)

if __name__ == '__main__':
    main(sys.argv)
