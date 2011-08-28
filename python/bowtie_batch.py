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
        single_end = line[2]
        indices_names = line[3:]
        indices = []
        for index_name in indices_names:
            index_name = index_name.split("-")
            indices.append((index_name[0], index_name[1]))
        for index in indices:
            line_args = ['-d', date, '-s', sample, '-i', index]
            bowtie.bowtie(date, sample, single_end, index)

if __name__ == "__main__":
    main(sys.argv)
