#! /usr/bin/env python

import sys, os
import argparse
import bowtie

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', required=True, dest='manifest', help='list of files to map')
    args = parser.parse_args()
    manifest = open(args.manifest)
    for line in manifest:
        line = line.split()
        date = line[0]
        sample = line[1]
        index = line[2]
        line_args = ['-d', date, '-s', sample, '-i', index]
        bowtie.bowtie(date, sample, index)

if __name__ == "__main__":
    main(sys.argv)
