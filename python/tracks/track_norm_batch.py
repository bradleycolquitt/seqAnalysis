#! /usr/bin/env python

import sys, os
import argparse

import genome.db

import track_norm as tn
from multiprocessing import Pool
from functools import partial

def worker(track_name, gdb=None, assembly=None, norm_type=None):
    tn.run(track_name, gdb, norm_type)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', dest='in_dir')
    parser.add_argument('-n', dest='norm_type', choices=['rpm'], default='rpm')
    parser.add_argument('-c', dest='ncores', type=int, help='number of threads to use', default=4)
    parser.add_argument('--assembly', help="assembly to use", default=None)
    args = parser.parse_args()

    gdb = genome.db.GenomeDB()
    all_files = gdb.list_tracks(subdir=args.in_dir, recursive=False)

    pool = Pool(processes=args.ncores)
    partial_worker = partial(worker,
                             gdb=gdb,
                             assembly=args.assembly,
                             norm_type=args.norm_type)

    partial_worker(all_files[0])
    #pool.map(partial_worker, all_files, 1)
    #pool.close()
    # for i in xrange(len(out_files)):
    #    print out_files[i]
    #    tracks = [ trim_files[i] + "-" + p for p in suffixes]
    #    combine_tracks.create_combined_tracks(out_files[i], tracks, args.assembly,
    #                       np.dtype(args.dtype))


if __name__ == '__main__':
    main()
