#!/usr/bin/env pytho

import sys
import os
import argparse
import stats_bf

FEATURE_PATH = '/home/user/lib/features_general'
OUT_PATH = '/media/storage2/analysis/features/bf'
SUMMARY_PATH_NORM = "/media/storage2/analysis/features/norm/summaries"
#UMMARY_PATH_UNNORM = 
SUMMARY_PATH_RAW = "/media/storage2/analysis/features/summaries"

def permute(l):
    one = 0
    ind = len(l)
    result = []
    while one < ind:
        two = one + 1
        while two < ind:
            result.append([l[one], l[two]])
            two = two + 1
        one = one + 1
    return result

def get_permute(path, file):
    sample_name = "/".join([path, file])
    sample_data = open(sample_name)
    samples = sample_data.readline().strip().split()
    sample_data.close()
    return(permute(samples))

def main(argv): 
    parser = argparse.ArgumentParser()
    parser.add_argument('--set')
    parser.add_argument('--feature')
    parser.add_argument('--feature_set', required=False)
    parser.add_argument('--data_type')
    parser.add_argument('--cutoff', type=float)
    args = parser.parse_args()
    
    path = ""
    if args.data_type == "norm":
        path = SUMMARY_PATH_NORM
    elif args.data_type == "unnorm":
        path = SUMMARY_PATH_NORM
    elif args.data_type == "raw":
        path = SUMMARY_PATH_RAW
    
    #if args.feature_set:
    #    files = os.listdir(path)
    #    files = filter(re.compile("args.set").search, files)
    #    for file in files:
    #        samples_permute = get_permute(path, file)
    #        ds = []
    #        for sample in samples_permute:
    #            ds.append(stats_bf.import_features(args.set, args.feature, sample, args.data_type))
    #        run(ds, out_path, samples_permute)
            
    sample_name = "/".join([path, "_".join([args.set, args.feature])])
    sample_data = open(sample_name)
    samples = sample_data.readline().strip().split()
    samples_permute = permute(samples)
    sample_data.close()
    print samples_permute
    ds = []
    for sample in samples_permute:
        ds.append(stats_bf.import_features(args.set, args.feature, sample, args.data_type))
        
    for i in range(len(ds)):
        n = ds[i].keys()
        bf = stats_bf.worker(ds[i])  
        out_path = "/".join([OUT_PATH, "_".join([args.set, args.feature, args.data_type])])
        if not os.path.exists(out_path): os.makedirs(out_path)
        out_name = "/".join([out_path, "_".join(samples_permute[i])])
        outfile = open(out_name, 'w')
        for i in range(len(bf)):
            if args.cutoff:
                if bf[i] > args.cutoff or bf[i] < -1 * args.cutoff:
                    outfile.write("%s\t%d\n" % (n[i], bf[i]))
            else:  
                outfile.write("%s\t%d\n" % (n[i], bf[i]))
        outfile.close()
         
if __name__ == '__main__':
    main(sys.argv)