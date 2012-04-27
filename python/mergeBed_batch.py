#!/usr/bin/env python

import sys, os, re
from subprocess import *

feature_path = "/home/user/lib/features_trim"

def main(argv):
    bed_dir = argv[1]
    beds = os.listdir(bed_dir)
    #feature_info_path = "/home/user/lib/feature_info"
    #if not os.path.exists(feature_info_path): os.makedirs(feature_info_path)
    
    #feature_span_out = open("/".join([feature_info_path, "feature_spans"]), 'w')
    #out = "feature\tfeature_span\ttotal_span\tfraction\n"
    #feature_span_out.write(out)
    
    #feature_spans = {}   # Set up dict of feature spans
    
    for bed in beds:
        if not re.search("merged", bed):
            print bed
            bed_file = "/".join([bed_dir, bed])
            merged_file = "_".join([bed_file, "merged"])
            merged_fileno = open(merged_file, "w")
            cmd_args = ['mergeBed', '-nms', '-i', bed_file]
            p = Popen(cmd_args, stdout=merged_fileno)
            p.wait()
        #feature_spans[feature] = 0
        #feature_file = open("/".join([feature_path, feature]))
        #for line in feature_file:
        #    line = line.split()
        ##    span = int(line[2]) - int(line[1])
        #   feature_spans[feature] = feature_spans[feature] + span
       
    ##span_total = sum(feature_spans.values())
    #for item in feature_spans.iteritems():
    #    fraction = float(item[1]) / span_total
    #    out = "{0}\t{1}\t{2}\t{3}\n".format(item[0], str(item[1]), str(span_total), str(fraction))
    #    feature_span_out.write(out)
    
if __name__ == '__main__':
    main(sys.argv)
