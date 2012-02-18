#! /usr/bin/env python

import sys, os
import argparse
import tophat

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
        errorlog = open("tophat_log", 'w')
        for index in indices:
            print line
            #line_args = ['-d', date, '-s', sample, '-i', index]
            if single_end == "PE": single_end = ""
            try:
                tophat.tophat(date, sample, bool(single_end), index)
            except:
                errorlog.write(str(sys.exc_info()[0]) + " ".join([str(index[0]), index[1]]) + "\n")
                continue
        errorlog.close()   

if __name__ == "__main__":
    main(sys.argv)
