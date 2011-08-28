#!/usr/bin/env python

import sys
import argparse
import tables as tb
from subprocess import Popen

TRACK_DIR = "/media/storage2/data/h5"
WIG_DIR = "/media/storage2/data/wig"

def main(argv):

    parser = argparse.ArgumentParser()
    parser.add_argument('-t', dest='track_file')
    args = parser.parse_args()
    track_file = TRACK_DIR + "/" + args.track_file
    h5 = tb.openFile(track_file, 'r')
    nodes = h5.iterNodes("/")
    nodes_tbp = []
    for node in nodes:
        node_name = node._v_name
        cmd_args = ['TrackToWig', '-t', track_file, '-n', node_name, '-w', WIG_DIR + "/" + node_name]
        print cmd_args
        p = Popen(cmd_args)
        p.wait()
        
if __name__ == "__main__":
    main(sys.argv)