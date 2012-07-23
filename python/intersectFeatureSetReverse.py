#!/usr/bin/env python

import sys, os
import argparse
from subprocess import *

feature_path = "/home/user/lib/features_merged"
out_path = "/media/storage2/data/homer/peaks/intersections_reverse"

def main(argv):
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', dest='reference_bed')
    args = parser.parse_args()
    
    #out_file = "/".join([out_path, "_".join([args.reference_bed, "features_general"])])
    features = os.listdir(feature_path)
    
    #Prepare output directory
    reference_out_path = "/".join([out_path, os.path.basename(args.reference_bed)])
    if not os.path.exists(reference_out_path): os.mkdir(reference_out_path)
    
    #Get total count
    #total_records = 0
    #cmd_args = ['wc', outfile]
    #p = Popen(cmd_args, stdout=total_records)
    #p.wait()
    #print total_records
    
    # Get length of reference bed
    cmd_args = ['wc', args.reference_bed]
    reference_length = Popen(cmd_args, stdout=PIPE).communicate()[0]
    reference_length = int(reference_length.split()[0])
    
    #p.wait()
    #print reference_length
    #return
        
    feature_length = 0
    count = 0
    
    # Set up summary file
    summary_file = open("/".join([reference_out_path, "summary"]), 'w')
    out = "feature\treference_length\tfeature_length\tcount\n"
    summary_file.write(out)
    
    # run intersectBed for reference_bed against each feature in features_general
    for feature in features:
        print feature
        feature_file_path = "/".join([feature_path, feature])
        feature_out_file = open("/".join([reference_out_path, feature]), 'w')
        cmd_args = ['intersectBed', '-b', args.reference_bed, '-a', feature_file_path]
        p = Popen(cmd_args, stdout=feature_out_file)
        p.wait()
       
        
    
    
if __name__ == '__main__':
    main(sys.argv)
