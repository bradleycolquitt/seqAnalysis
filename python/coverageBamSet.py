#!/usr/bin/env python

import sys, os, re
import argparse
import pysam
import pdb
from subprocess import *

feature_path = "/home/user/lib/features_merged"
out_path = "/media/storage2/analysis/features/coverage"



def main(argv):
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-bam', dest='bam')
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
    bam_prefix = os.path.basename(args.bam).split(".bam")[0]
    reference_out_path = "/".join([out_path, bam_prefix])
    if not os.path.exists(reference_out_path): os.mkdir(reference_out_path)
    
    #Get total count
    #total_records = 0
    #cmd_args = ['wc', outfile]
    #p = Popen(cmd_args, stdout=total_records)
    #p.wait()
    #print total_records
    
    # Get total reads
    bamfile = pysam.Samfile(args.bam, 'rb')
    nreads = bamfile.mapped
    
    # Get length of reference bed
    #cmd_args = ['wc', args.reference_bed]
    #reference_length = Popen(cmd_args, stdout=PIPE).communicate()[0]
    #reference_length = int(reference_length.split()[0])
    
    #p.wait()
    #print reference_length
    #return
        
    #feature_length = 0
    #span = 0
    #count = 0
    #pdb.set_trace()
    # Get average length of peaks
    #peak_lengths = np.zeros(reference_length)
    #peak_lengths_cum = 0
    #reference_file = open(args.reference_bed)
    #for line in reference_file:
    #    line = line.split()
    #    peak_lengths_cum = peak_lengths_cum + (int(line[2]) - int(line[1]))
    #peak_lengths_mean = peak_lengths_cum / reference_length 
    
    # Set up summary file
    #summary_file = open("/".join([reference_out_path, "summary"]), 'w')
    #out = "feature\treference_length\tfeature_length\tfeature_span\tcount\tfraction\tfraction_norm\n"
    #summary_file.write(out)
    
    # run intersectBed for reference_bed against each feature in features_general
    for feature in features:
        print feature
        #pdb.set_trace()
        feature_file_path = "/".join([feature_path, feature])
        feature_out_file_name = "/".join([reference_out_path, feature])
        feature_out_file = open(feature_out_file_name, 'w')
        cmd_args = ['multiBamCov', '-D', '-bams', args.bam, '-bed', feature_file_path]
        p = Popen(cmd_args, stdout=feature_out_file)
        p.wait()
        feature_out_file.close()
        #feature_out_file = open("/".join([reference_out_path, feature]), 'r')
        
        # Get length feature bed
        #cmd_args = ['wc', feature_file_path]
        #feature_length = Popen(cmd_args, stdout=PIPE).communicate()[0]
        #feature_length = int(feature_length.split()[0])
        
        # Get span feature bed
        feature_file = open(feature_file_path)
        for line in feature_file:
            line = line.split()
            span = span + (int(line[2]) - int(line[1])) #+ (2 * (peak_lengths_mean - 1))        #pdb.set_trace()
            
        
        #reference_[feature] = 0
        #feature_length[feature] = 0
        feature_out_file_mod_name = "/".join([reference_out_path, feature]) + "mod"
        feature_out_file_mod = open(feature_out_file_mod_name, 'w')
        for line in feature_out_file:
            line = line.split()
            span = line[2] - line[1]
            count = line[6]
            count = (1E6 * count) / (span * nreads)
            out = "/t".join([line[:6], count] + "\n")
            feature_out_file_mod.write(out)
            #total_overlap[feature] = total_overlap[feature] + int(line[5])
            #total_length[feature] = total_length[feature] + int(line[6])
            #feature_length = feature_length + 1
            #count = count + int(line[6])
            #fraction = round(float(total_overlap[feature]) / float(total_length[feature]), 3)
        #fraction = float(count) / reference_length
        #fraction_norm = 1000000 * fraction / span
            #out = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(feature, str(reference_length), str(feature_length),
            #                                               str(span), str(count), str(fraction), str(fraction_norm))
#        summary_file.write(out)
#        count = 0
#        span = 0
        feature_out_file.close()
        feature_out_file_mod.close()
        os.rename(feature_out_file_mod_name, feature_out_file_name)
        
    
    
if __name__ == '__main__':
    main(sys.argv)
