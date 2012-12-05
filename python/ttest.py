#!/usr/bin/env python

import sys
import os

FEATURE_PATH = "/media/storage2/analysis/features/norm"

def import_data(feature, samples):
    sample_names = [FEATURE_PATH]
    for s1, s2 in itertools.izip(sample_data[0], sample_data[1]):
        s1 = s1.split()
        s2 = s2.split()
        d[s1[3]] = (float(s1[4]), float(s2[4]))


def main(argv):
    data = import_data()
    
        
if __name__ == '__main__':
    main(sys.argv)