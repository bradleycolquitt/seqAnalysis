#!/usr/bin/env python


import sys
import os
import re
import pdb
import argparse
import numpy as np

import genome.db
import combine_tracks

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', dest='in_dir')
    parser.add_argument('-p', nargs='?', dest='suffixes', help="suffixes of files to combine")
    parser.add_argument("--dtype", metavar="", action="store",
                        choices=("uint8", "uint16"), default="uint8",
                        help="datatype of combined track")

    parser.add_argument('--assembly', help="assembly to use", default=None)
    args = parser.parse_args()

    gdb = genome.db.GenomeDB()
    all_files = gdb.list_tracks(subdir=args.in_dir, recursive=False)

    if args.suffixes is None:
        suffixes = ['fwd', 'rev']
    pattern = "-(" + "|".join(suffixes) + ")"

    trim_files = list(set([ re.sub(pattern, '', f) for f in all_files ]))
    base_files = [os.path.basename(f) for f in trim_files]

    out_files = ["/".join([args.in_dir, 'combined', f]) for f in base_files]
    #pdb.set_trace()

    for i in xrange(len(out_files)):
        print out_files[i]
        tracks = [ trim_files[i] + "-" + p for p in suffixes]
        combine_tracks.create_combined_tracks(out_files[i], tracks, args.assembly,
                           np.dtype(args.dtype))


if __name__ == '__main__':
    main()
