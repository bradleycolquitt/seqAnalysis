#! /usr/bin/env python

#### 
# Find position of through between signal peaks
# INPUT: signal matrix (rows are features, columns are positions)
# INPUT: BED file with column 4 matching column 1 of signal matrix
# OUTPUT: BED file with start and end adjusted so that midpoint is trough; length is the same as input
#
# General strategy:
#     Cross-correlate input signal with template of two gaussian
#     Vary intergaussian distance, collect maximum correlation value and position
#     Loop through input lines
#     Loop through range of range of bandwiths for 
####

import sys
import pdb
import itertools
import numpy as np
import scipy.signal as ss

def double_gaussian(M, d, j):
    #pdb.set_trace()
    vect = np.zeros(M)
    vect_mid = M / 2
    gaus = ss.gaussian(j, std=1)

    start = vect_mid - ((d / 2) + (j/2))
    end = vect_mid - ((d / 2) - (j/2))
    vect[start:end] = vect[start:end] + gaus

    start = vect_mid + ((d / 2) - (j/2))
    end = vect_mid + ((d / 2) + (j/2))
    vect[start:end] = vect[start:end] + gaus

    return vect

def cross_gauss(input, template):    
    out = ss.correlate(input, template)

def adjust_bed(bed_row, adjust):
    bed_row = bed_row.strip().split()
    start = int(bed_row[1])
    end = int(bed_row[2])
    mid = (start + end) / 2
    adjust = mid - adjust
    bed_row[1] = str(start - adjust)
    bed_row[2] = str(end - adjust)
    bed_row = "\t".join(bed_row) + "\n"
    return bed_row

def main(argv):

    pdb.set_trace()
    # Load input matrix and bed
    input = open(argv[1])
    bed = open(argv[2])
    outfile = open(argv[3], 'w')
    
    # Define range of peak distances
    peak_distances = xrange(5, int(argv[4]))

    # Loop through input matrix and bed
    for (input_row, bed_row) in itertools.izip(input, bed):
        print(input_row.split()[0], bed_row.split()[3])
        assert(input_row.split()[0] == bed_row.split()[3])
        input_vect = np.array(input_row.strip().split()[2:], dtype="float")
        distance_values = []
        # Loop through peak distances
        for distance in peak_distances:
            # Create double gaussian
            input_length = len(input_row)
            gaus = double_gaussian(input_length, distance, 10)
            
            # Correlate with input vector
            corr_vector = ss.correlate(input_vect, gaus)
            
            # Add maximum value and position (as tuple) to list
            max_value = max(corr_vector)
            max_pos = [i for i, j in enumerate(corr_vector) if j == max_value]
            distance_values.append((max_value, max_pos[0]))
   
        # Find max correlation and position
        max_value = (0,0)
        for i,j in enumerate(distance_values):
            if j > max_value[1]:
                max_value = (i,j)

        # Adjust bed to center found position
        out_bed = adjust_bed(bed_row, max_value[0])
        outfile.write(out_bed)
    

if __name__ == "__main__":
    main(sys.argv)
