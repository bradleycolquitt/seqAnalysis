#!/usr/bin/env python

import sys, os
from subprocess import *

feature_path = "/home/user/lib/features_merged"

def main(argv):
    
    features = os.listdir(feature_path)
    
    feature_info_path = "/home/user/lib/feature_info"
    if not os.path.exists(feature_info_path): os.makedirs(feature_info_path)
    
    feature_span_out = open("/".join([feature_info_path, "feature_merged_spans"]), 'w')
    out = "feature\tfeature_span\ttotal_span\tfraction\n"
    feature_span_out.write(out)
    
    feature_spans = {}   # Set up dict of feature spans
    
    for feature in features:
        print feature
        feature_spans[feature] = 0
        feature_file = open("/".join([feature_path, feature]))
        for line in feature_file:
            line = line.split()
            span = int(line[2]) - int(line[1])
            feature_spans[feature] = feature_spans[feature] + span
       
    span_total = sum(feature_spans.values())
    for item in feature_spans.iteritems():
        fraction = float(item[1]) / span_total
        out = "{0}\t{1}\t{2}\t{3}\n".format(item[0], str(item[1]), str(span_total), str(fraction))
        feature_span_out.write(out)
    
if __name__ == '__main__':
    main(sys.argv)
