#! /usr/bin/env python

import sys, os, re
import argparse
import tophat2

# Process manifest file for tophat2 mapping
# Format: date sample SE|PE mean_size sd_size GTF library-type species

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', required=True, dest='manifest', help='list of files to map')
    args = parser.parse_args()
    manifest = open(args.manifest)
    for line in manifest:
        if not re.search("#", line):
            line = line.split()
            date = line[0]
            sample = line[1]
            single_end = line[2]
            mean = line[3]
            sd = line[4]
            gtf = line[5]
            library_type = line[6]
            species = line[7]
            #indices_names = line[8:]
            #indices = []
            #for index_name in indices_names:
            #    print index_name
            #    index_name = index_name.split("-")
            #    indices.append((index_name[0], index_name[1]))
            errorlog = open("tophat_log", 'w')
            #for index in indices:
                #print index
                #line_args = ['-d', date, '-s', sample, '-i', index]
            if single_end == "PE": single_end = ""
            try:
                tophat2.tophat(date, sample, bool(single_end), mean, sd, gtf, library_type, species)
            except:
                errorlog.write(str(sys.exc_info()[0]) + sample + "\n")
                raise
                    
            errorlog.close()   

if __name__ == "__main__":
    main(sys.argv)
