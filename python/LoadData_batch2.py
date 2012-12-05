#!/usr/bin/env python

import sys
import argparse
import tables as tb
from subprocess import Popen

TRACK_DIR = "/media/storage2/data/h5"
WIG_DIR = "/media/storage2/data/wig"

## Manifest format
##  Field 1: file location
##  Field 2: track name

def main(argv):

    parser = argparse.ArgumentParser()
    parser.add_argument('-t', dest='track_file')
    parser.add_argument('-m', dest='manifest')
    parser.add_argument('-r', dest='resolution')
    args = parser.parse_args()
    track_file = TRACK_DIR + "/" + args.track_file
    manifest = open(args.manifest)
#    h5 = tb.openFile(track_file, 'w')
#    nodes = h5.iterNodes("/")
#    nodes_tbp = []
    for line in manifest:
        print line
        line = line.split()
        cmd_args = ['LoadData',
                    '-i', line[0],
                    '-o', track_file,
                    '-t', line[1],
                    '-g', '/media/storage2/genomedata/chromosomes.trk',
                    '-n', args.resolution]
        p = Popen(cmd_args)
        p.wait()
        
if __name__ == "__main__":
    main(sys.argv)