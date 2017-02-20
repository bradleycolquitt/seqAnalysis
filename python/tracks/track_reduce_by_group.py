#! /usr/bin/env python

import os
import sys
import argparse
import tables as tb
import numpy as np
import re
import pdb
#import track_util as tutil

import genome.track
import genome.db
import combine_tracks2 as ct

import pandas as pd

def read_df(df_fname):
    info = pd.DataFrame.from_csv(df_fname, sep='\t')

    return info

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', dest='track_dir', help='directory containing tracks')
    parser.add_argument('-i', dest='info', help='TSV containing metadata')
    parser.add_argument('-o', dest='out_name', help='out track file')
    parser.add_argument('-f', dest='function', default="mean")
    parser.add_argument('-g', dest='groups', default=None, required=True, nargs='*')

    parser.add_argument("--dtype", metavar="", action="store",
                        choices=("uint8", "uint16"), default="uint8",
                        help="datatype of combined track")

    parser.add_argument('--assembly', help="assembly to use", default=None)

    args = parser.parse_args()
    groups = args.groups

    # open genome
    gdb = genome.db.GenomeDB()
    tracks = gdb.list_tracks(subdir=args.track_dir, recursive=False)
    tracks_base = [os.path.basename(t) for t in tracks]
    out_prefix = "/".join([os.path.dirname(args.track_dir), "reduced", "-".join(groups)])

    # read info file and group
    info = read_df(args.info)
    info_gp = info.groupby(groups)

    for name, group in info_gp:
        # match filenames with groups
        curr_tracks = group['group_name']
        curr_tracks = [f.strip('/') for f in curr_tracks.tolist()]
        ind = [i for i,v in enumerate(tracks_base) if v in curr_tracks]
        curr_tracks_full = [tracks[i] for i in ind]

        if len(curr_tracks_full) == 0:
            continue

        if len(name) == 1:
            curr_name = name
        else:
            curr_name = "-".join(list(name))
        print curr_name
        out_fname = "/".join([out_prefix, curr_name])

        if gdb.has_track(out_fname):
            gdb.delete_track(out_fname)


        ct.create_combined_tracks(out_fname, curr_tracks_full, args.assembly, args.function,
                                  dtype=np.dtype(args.dtype))

if __name__ == '__main__':
    main(sys.argv)
