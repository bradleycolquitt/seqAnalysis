#!/usr/bin/env python

import sys
import argparse
import tables as tb
import track_util as tutil

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', dest='track', help='input track')
    parser.add_argument('-n', dest='subtrack', help="subtrack name. Input 'all' to process all subtracks")

    args = parser.parse_args()
    track_file = tb.openFile(args.track, 'a')
    subtrack_name = args.subtrack

    track_file.remove_node("/" + subtrack_name, recursive=True)
    track_file.flush()
    track_file.close()
if __name__ == '__main__':
    main(sys.argv)
