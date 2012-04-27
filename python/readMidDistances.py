#!/usr/bin/env python

####################################################################
# Compute variance of distances of the midpoint of reads in a given BAM #
# that overlap each genomic position                               #
#   Param: BAM file                                                #
#   Output: Track file with with (max(variance) - variance)/max(variance)     #
####################################################################

import sys, os
import argparse
import pdb
import pysam
import random
import copy
import numpy as np
from subprocess import Popen
from subprocess import PIPE

TEST_HEADER = { 'HD': {'VN': '1.0'},
                'SQ': [{'LN': 10000, 'SN': 'chr1'}] }

def construct_read(count, ref, pos, isize, is_read1):
    a = pysam.AlignedRead()
    a.qname = "read" + str(count)
    a.seq = "N"*50
    flag = 93
    if is_read1: flag = 63
    a.flag = flag
    a.tid = ref
    a.pos = pos
    a.mapq = 20
    a.is_duplicate = False
    a.is_paired = True
    a.is_proper_pair = True
    a.is_read1 = is_read1
    a.is_read2 = not is_read1
    a.is_reverse = not is_read1
    a.is_unmapped = False
    a.mate_is_reverse = is_read1
    a.mate_is_unmapped = False
    a.cigar = ((0,50),)
    a.rnext = 0
    a.pnext = pos + isize + 1
    a.mpos = pos + isize
    a.qual = "I"*50
    a.isize = isize
    a.tlen = isize
    return a

def make_test_bam(total_reads, jigg, out):
    # Start with initial position at chr1 200
    # First define positions by midpoint
    # That is subtract isize from defined initial position
    # Define isizes by normal distribution with mean 145, sd 10
    # Step by 200
    # Continue for 10000
    # For each midpoint make Read1, Read2 AlignedRead
    #   Negative isize for Read2
    test_bam = pysam.Samfile(out + ".bam", 'wb', header=TEST_HEADER)
    #pdb.set_trace()
    mid_position = 200
    #Randomly generate number from 2 to 20 for N value
    read_range = xrange(2, 20)
    count = 0
    while total_reads > 0:
        N = random.sample(read_range, 1)[0]
        total_reads = total_reads - N
        while N > 0:
            count = count + 1
            N = N - 1
            # Introduce jiggle
            tmp_mid_position = int(round(mid_position + random.uniform(-1, 1,) * jigg))
            isize = int(round(random.gauss(145, 10)))
            positions = ((tmp_mid_position - (isize / 2), isize), (tmp_mid_position + (isize / 2), -isize))
            #abs might not be right
            for position in positions:
            # Construct/write Aligned read
                read = construct_read(count, 0, position[0], position[1], position[1] > 0)
                test_bam.write(read)
        mid_position = mid_position + 300
    test_bam.close()
    #pysam.sort(out + ".bam", out + "_sort")
    #pysam.index(out+ "_sort.bam")
    pysam.index(out + ".bam")

def extend_bam(bam):
    bam_prefix = bam.split(".bam")[0]
    # bedtools bamToBed to convert to Bed

    #bed_name = bam_prefix + ".bed"
    #bed_file = open(bed_name, 'w')
    bam_frag_name = bam_prefix + "_frag.bam"
    bam_frag_file = open(bam_frag_name, 'w')
    cmd_args1 = ['bamToBed', '-i', bam, '-bedpe']
    cmd_args2 = ['cut', '-f', '1-2,6-9']
    cmd_args3 = ['bedToBam', '-i', 'stdin', '-g', '/seq/lib/mouse.mm9.genome_norand']
    p1 = Popen(cmd_args1, stdout=PIPE)
    p2 = Popen(cmd_args2, stdin=p1.stdout, stdout=PIPE)
    p3 = Popen(cmd_args3, stdin=p2.stdout, stdout=bam_frag_file)
    p3.wait()
    pysam.sort(bam_frag_name, bam_prefix + "_frag_sort")
    pysam.index(bam_prefix + "_frag_sort.bam")
    #bed_file.close()

    # cut columns 1,2,6,7,8,9 to make fragment bed
    #frag_name = bam_prefix + "_frag.bed"
    #frag_file = open(frag_name, 'w')
    
    #p = Popen(cmd_args, stdout=frag_file)
    #p.wait()
    #frag_file.close()
    #finally:
    #    os.remove(bed_name)
    #try:
        # convert fragment to bed to bam using bedToBam
        
        
        
    #    p = Popen(cmd_args, stdout=bam_frag_file)
    #    p.wait()
    #    bam_frag_file.close()
    #finally:
    #    os.remove(bam_frag_name)
    
    
def extract_worker(sam_file, h5_file, ref, ref_length, track_name):
    pass
# Take paired end BAM file
# For chr in reference list:
# Initialize numpy vector with length equal to chr length
# Generate pileup iterator for chr
# For pileupProxy in iterator:
#   For pileupRead.alignment in pileupProxy:
#       add absolute value of insert size to list of insert sizes
#   average sizes
#   write to value vector
# Write out as track in HDF5 file:

def main(argv):
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', dest='bam')
    parser.add_argument('-w', dest='wig')
    parser.add_argument('-r', dest='resolution')
    args = parser.parse_args()
    
    sam_file = pysam.Samfile(args.bam, 'rb')
    #h5_file = tb.openFile(output, 'a')
    wig_file = open(args.wig, 'w')
    refs = sam_file.references
    ref_lengths = sam_file.lengths
    resolution = int(args.resolution)
    #pool = Pool(processes=6)
    
    #start_time = time.time()
    #args = []
    #args = [(sam_file, h5_file, ref, ref_lengths[])for files in files_group]
    #for arg in args:
    #    pool.apply_async(extract_worker, (arg,))
    #    index_split(files, argv)
    #pool.close()
    #pool.join()
    ## Loop through chromosomes
    for chr_index in range(len(refs)):
        ## Extract and initialize variables
        curr_ref = refs[chr_index]
        print curr_ref
        #pdb.set_trace()
        curr_length = ref_lengths[chr_index]
        if resolution > 1:
            curr_length = curr_length / resolution
        values = np.zeros(curr_length)
        distances = np.zeros(1000)
        distances_array_len = len(distances)
        res_count = 0
        insert_size = 0
        sum_distances = 0
        num_reads = 0
        #value_ind = 0
        position = 0
        position_mean = 0
        position_var = 0
        sum_position_mean = 0
        sum_position_var = 0
        total_num_reads = 0
        value_ind = 0
        value_out = 0
        first_pos = True
        #read_ind = 0
        
        out = "fixedStep chrom={0} start=1 step={1} span={1}\n".format(curr_ref, resolution)
        wig_file.write(out)
        
        ## Create pileup iterator
        it = sam_file.pileup(reference=curr_ref)
        
        ## Loop through each position in iterator
        print "Computing variances"
        for proxy in it:
            position = proxy.pos
            #if first_pos:
            value_ind = position / resolution
            #    first_pos = False
            #if res_count == 0: value_ind = value_ind + 1
            #assert(value_ind < curr_length)
            
            pdb.set_trace()
            ## Loop through each read at current position
            for pread in proxy.pileups:
                insert_size = pread.alignment.isize
                if insert_size > 50:
                    mid_distance = abs((pread.alignment.pos + (insert_size/2)) - position)
                    distances[num_reads] = mid_distance
                    #distances.append(mid_distance)
                    #sum_distances = sum_distances + abs(read_mid)
                    num_reads = num_reads + 1
                    if num_reads == distances_array_len: break
                
            if num_reads > 5:
                #pdb.set_trace()
                position_mean = np.mean(distances[:num_reads])
                position_var = sum(pow((distances[:num_reads] - position_mean), 2)) / (num_reads - 1)
                total_num_reads = total_num_reads + num_reads
            sum_position_mean = sum_position_mean + position_mean
            sum_position_var = sum_position_var + position_var
            #pdb.set_trace()
            #values[position] = position_var
            #print position_mean, position_var
            #wig_file.write(str(value) + "\n")
            #distances = []
            num_reads = 0
            #value_ind = value_ind + 1
            position_mean = 0
            position_var = 0
            res_count = res_count + 1   
            if res_count == resolution:
                #pdb.set_trace()
                if sum_position_var > 0:
                    #value_out = (sum_position_mean / sum_position_var) / resolution
                    #value_out = 1 / ((1 + sum_position_mean) * sum_position_var * resolution) 
                    value_out = total_num_reads / (sum_position_var * resolution)
                    #value_out = sum_position_mean / resolution
                    #value_out = sum_position_var / resolution
                values[value_ind] = value_out
                #value_ind = value_ind + 1
                value_out = 0
                res_count = 0
                sum_position_mean = 0
                sum_position_var = 0
                total_num_reads = 0       
        ## Write HDF5 track
        #print "Writing"
        #wig_file.write()
        #track_util.writeTrack(h5_file, track_name, curr_ref, values, 1)
        #print "Computing scores"
        max_value = np.max(values)
        #print "Normalizing"
        scores = [0]
        
        if max_value > 0:
            scores = values
        #    scores = (max_value - values) / max_value
        
        print "Writing scores"
        #for value in values:
        #    print>>wig_file, str(value)
        for score in scores:
            print>>wig_file, str(score)
            #wig_file.write(str(score) + "\n")
            
    #cmd_args = ['igvtools', 'tile', output, output + ".tdf", 'mm9']
    ##p = Popen(cmd_args)
    #p.wait()
    #h5_file.flush()
    #h5_file.close()
if __name__ == '__main__':
    main(sys.argv)