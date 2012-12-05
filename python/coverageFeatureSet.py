#!/usr/bin/env python

import sys, os
import argparse
from subprocess import Popen

feature_path = "/home/user/lib/features_trim"
out_path = "/media/storage2/data/homer/peaks/coverage"

def main(argv):
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', dest='reference_bed')
    args = parser.parse_args()
    
    #out_file = "/".join([out_path, "_".join([args.reference_bed, "features_general"])])
    features = os.listdir(feature_path)
    
    #Prepare output directory
    reference_out_path = "/".join([out_path, args.reference_bed])
    if not os.path.exists(reference_out_path): os.mkdir(reference_out_path)
    
    #Get total count
    #total_records = 0
    #cmd_args = ['wc', outfile]
    #p = Popen(cmd_args, stdout=total_records)
    #p.wait()
    #print total_records
    
    # Define dicts for total coverage and total length for each feature
    total_overlap = {}
    total_length = {}
    
    # run coverageBed for reference_bed against each feature in features_general
    
    summary_file = open("/".join([reference_out_path, "summary"]), 'w')
    out = "feature\toverlap\tlength\tfraction\n"
    summary_file.write(out)
    
    for feature in features:
        print feature
        feature_out_file = open("/".join([reference_out_path, feature]), 'w')
        cmd_args = ['coverageBed', '-a', args.reference_bed, '-b', "/".join([feature_path, feature])]
        p = Popen(cmd_args, stdout=feature_out_file)
        p.wait()
        feature_out_file.close()
        feature_out_file = open("/".join([reference_out_path, feature]), 'r')

        total_overlap[feature] = 0
        total_length[feature] = 0
        for line in feature_out_file:
            line = line.split()
            total_overlap[feature] = total_overlap[feature] + int(line[5])
            total_length[feature] = total_length[feature] + int(line[6])
        fraction = round(float(total_overlap[feature]) / float(total_length[feature]), 3)
        out = "{0}\t{1}\t{2}\t{3}\n".format(feature, total_overlap[feature], total_length[feature], fraction)
        summary_file.write(out)
        
        
    
    
if __name__ == '__main__':
    main(sys.argv)
