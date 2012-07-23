#!/usr/bin/env python

import sys, os, re
import argparse
from subprocess import *

feature_path = "/home/user/lib/features_merged"
out_path = "/media/storage2/data/homer/peaks/intersections"



def main(argv):
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', dest='reference_bed')
    parser.add_argument('--type', dest='feature_type')
    #parser.add_argument('-d', dest='reference_dir')
    args = parser.parse_args()
    
    if args.feature_type == "rmsk":
        feature_path = "/home/user/lib/features_rmsk_merged"
    else:
        feature_path = "/home/user/lib/features_merged"
    #out_file = "/".join([out_path, "_".join([args.reference_bed, "features_general"])])
    features = os.listdir(feature_path)
    features = [f for f in features if not re.search("not_used", f)]
    
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
    span = 0
    count = 0
    #pdb.set_trace()
    # Get average length of peaks
    #peak_lengths = np.zeros(reference_length)
    peak_lengths_cum = 0
    reference_file = open(args.reference_bed)
    for line in reference_file:
        line = line.split()
        peak_lengths_cum = peak_lengths_cum + (int(line[2]) - int(line[1]))
    peak_lengths_mean = peak_lengths_cum / reference_length 
    
    # Set up summary file
    summary_file = open("/".join([reference_out_path, "summary"]), 'w')
    out = "feature\treference_length\tfeature_length\tfeature_span\tcount\tfraction\tfraction_norm\n"
    summary_file.write(out)
    
    # run intersectBed for reference_bed against each feature in features_general
    for feature in features:
        print feature
        feature_file_path = "/".join([feature_path, feature])
        feature_out_file = open("/".join([reference_out_path, feature]), 'w')
        cmd_args = ['intersectBed', '-a', args.reference_bed, '-b', feature_file_path, '-c',
                    '-f', '0.5']
        p = Popen(cmd_args, stdout=feature_out_file)
        p.wait()
        feature_out_file.close()
        feature_out_file = open("/".join([reference_out_path, feature]), 'r')
        
        # Get length feature bed
        cmd_args = ['wc', feature_file_path]
        feature_length = Popen(cmd_args, stdout=PIPE).communicate()[0]
        feature_length = int(feature_length.split()[0])
        
        # Get span feature bed
        feature_file = open(feature_file_path)
        for line in feature_file:
            line = line.split()
            span = span + (int(line[2]) - int(line[1])) #+ (2 * (peak_lengths_mean - 1))        #pdb.set_trace()
            
        
        #reference_[feature] = 0
        #feature_length[feature] = 0
        for line in feature_out_file:
            line = line.split()
            #total_overlap[feature] = total_overlap[feature] + int(line[5])
            #total_length[feature] = total_length[feature] + int(line[6])
            #feature_length = feature_length + 1
            count = count + int(line[6])
            #fraction = round(float(total_overlap[feature]) / float(total_length[feature]), 3)
        fraction = float(count) / reference_length
        fraction_norm = 1000000 * fraction / span
        out = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(feature, str(reference_length), str(feature_length),
                                                           str(span), str(count), str(fraction), str(fraction_norm))
        summary_file.write(out)
        count = 0
        span = 0
        
        
    
    
if __name__ == '__main__':
    main(sys.argv)
