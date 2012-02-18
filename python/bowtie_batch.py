#! /usr/bin/env python

import sys, os
import argparse
import bowtie
import datetime

## Take manifest file of fastq files to map with bowtie
## Each line in manifest file has the form:
##      date [111221]   lane [s_1]  library-type [SE | PE] index-sample_name [1-moe_native_pol2]

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', required=True, dest='manifest', help='list of files to map')
    args = parser.parse_args()
    manifest = open(args.manifest)
    now = datetime.datetime.now()
    errorlog = open("bowtie_log" + "_" + str(now.year) + str(now.month) + str(now.day) + "_" + str(now.hour) + str(now.minute), 'w')
    for line in manifest:
        line = line.split()
        date = line[0]  ## date of sequencing
        lane = line[1]  ## lane of flow cell
        single_end = line[2] ## single or paired end
        indices_names = line[3:] ## list of indices to process
        indices = []
        for index_name in indices_names:
            index_name = index_name.split("-")
            indices.append((index_name[0], index_name[1]))
       
        
        for index in indices:
            print line
            if single_end == "PE": single_end = ""
            try:
                bowtie.bowtie(date, lane, bool(single_end), index)
            except:
                errorlog.write(str(sys.exc_info()[0]) + " ".join([str(index[0]), index[1]]) + "\n")
                continue
        errorlog.close()   

if __name__ == "__main__":
    main(sys.argv)
