#! /usr/bin/env python

import sys, os, re
import argparse
import bowtie2
import datetime

## Take manifest file of fastq files to map with bowtie
## Each line in manifest file has the form:
##      date [111221]   sample_directory  library-type [SE | PE] index-sample_name [1-moe_native_pol2]

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', required=True, dest='manifest', help='list of files to map')
    args = parser.parse_args()
    manifest = open(args.manifest)
    now = datetime.datetime.now()
    if not os.path.exists("logs"): os.mkdir("logs")
    errorlog = open("logs/bowtie_log" + "_" + str(now.year) + str(now.month) +
                    str(now.day), 'wa')
    
    for line in manifest:
        if not re.search("#", line):
            line = line.split()
            date = line[0]  ## date of sequencing
            sample = line[1] ## sample directory
            single_end = line[2] ## single or paired end
            #style = line[3]
            subsamples = line[3:] ## list of samples within directory to process
            #indices = []
            
            #for index_name in indices_names:
            #    index_name = index_name.split("-")
            #    indices.append((index_name[0], index_name[1]))
        
            for subsample in subsamples:
                print line
                if single_end == "PE": single_end = ""
                try:
                    bowtie2.bowtie(date, sample, bool(single_end), subsample)
                except:
                    errorlog.write(str(sys.exc_info()[0]) + " ".join([sample, subsample]) + "\n")
                    continue
    errorlog.close()   

if __name__ == "__main__":
    main(sys.argv)
