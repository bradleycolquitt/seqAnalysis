#!/usr/bin/env python

#########################################################################
# Compute variance of distances of the midpoint of reads in a given BAM #
# that overlap each genomic position                                    #
#   Param: BAM file                                                     #
#   Output: Wig file with with num_of_reads / variance                 #
########################################################################

import sys, os
import re
import shutil
import argparse
import pdb
import pysam
import random
import sequtils
import numpy as np
import multiprocessing
from subprocess import Popen
from subprocess import PIPE

TEST_HEADER = { 'HD': {'VN': '1.0'},
                'SQ': [{'LN': 500000, 'SN': 'chr1'}] }

SIM_NUC_PATH = "/home/user/data/nuc/sim"

class simulate_nuc:
    def __init__(self, total_reads, jitter, insert_sd, advance, jitter_advance, out):
        self.total_reads = total_reads
        self.jitter = jitter
        self.insert_sd = insert_sd
        self.advance = advance
        self.jitter_advance = jitter_advance
        self.out = out
        self.initial_mid_position = 200
        
    def run(self):
        sim_bam = pysam.Samfile(self.out + ".bam", 'wb', header=TEST_HEADER)
        #Randomly generate number from 2 to value normalized for genomic region for N value
        read_range = xrange(2, self.total_reads/500)
        #pdb.set_trace()
        mid_position = self.initial_mid_position

        count = 0
        total_reads = self.total_reads
        while total_reads > 0:
            N = random.sample(read_range, 1)[0]
            total_reads = total_reads - N
            while N > 0:
                count = count + 1
                N = N - 1
            
                # Introduce jitter
                tmp_mid_position = int(round(mid_position + random.uniform(-1, 1,) * self.jitter))
            
                isize = int(round(random.gauss(145, self.insert_sd)))
                positions = ((tmp_mid_position - (isize / 2), isize),
                             (tmp_mid_position + (isize / 2), -isize))
        
                # Construct/write paired AlignedReads           
                for position in positions:
                    read = construct_read(count, 0, position[0], position[1], position[1] > 0)
                    sim_bam.write(read)
                
                # Advance position by 300
            advance = self.advance
            if self.jitter_advance > 0: advance = int(round(random.gauss(advance, self.jitter_advance)))
            mid_position = mid_position + advance
        sim_bam.close()

class desc_file:
    def __init__(self, path):
        self.path = path
        self.file = open(self.path + "/desc", 'w')
    def write(self, total_reads, B, jitter_range, insert_size, jitter_advance):
        self.file.write("Simulated nucleosomes\n")
        self.file.write("Total number of reads: " + str(total_reads) + "\n")
        self.file.write("Number of simulations: " + str(B) + "\n")
        self.file.write("Range of midpoint jitterle factors: " + str(" ".join(str(jitter_range))) + "\n")
        self.file.write("Standard deviation of insert sizes (mean 145): " + str(insert_size) + "\n")
        if jitter_advance > 0:
            self.file.write("Nucleosome advanced by draw from normal distribution with mean 300, sd %s\n" % str(jitter_advance))

# Simulate B multiple nucleosome BAMs over given set of jitters (jitter_range)         
def simulate_nuc_range_multiple(jitter_range, set_no, total_reads=10000,
                                 B=100, insert_sd=10, advance=300, jitter_advance=0):

    read_path = "/".join([SIM_NUC_PATH, "set" + str(set_no)])
    if not os.path.exists(read_path): os.makedirs(read_path)
    desc = desc_file(read_path)
    desc.write(total_reads, B, jitter_range, insert_sd, jitter_advance)
    
    for jitter in jitter_range:
        print jitter
        jitter_path = "/".join([read_path, str(jitter)])
        if not os.path.exists(jitter_path): os.mkdir(jitter_path)
        pool = multiprocessing.Pool(processes=8)    
        args = []
        for b in xrange(B):
            out_bam = "/".join([jitter_path, "sim" + str(b)])
            args = [total_reads, jitter, insert_sd, advance, jitter_advance, out_bam]
            pool.apply_async(sim_measure_worker, (args,))
            #pool.apply(sim_measure_worker, (args,))
            #sim_measure_worker(args)
        pool.close()
        pool.join()
        
def sim_measure_worker(args):
    sim = simulate_nuc(args[0], args[1], args[2], args[3], args[4], args[5])
    #print "Sim run"
    sim.run()
    #print "Extend"
    extend_bam(args[5] + ".bam")
    #print "Distance"
    calc_distance_score(args[5] + "_frag_sort.bam", args[5] + "_frag", 1, False)
          
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
    bam_frag_file.close()
    pysam.sort(bam_frag_name, bam_prefix + "_frag_sort")
    
    cmd_args1 = ['samtools', 'view', '-h', bam]
    cmd_args2 = ['samtools', 'reheader', '-', bam_prefix + "_frag_sort.bam"]
    tmp = open(bam_prefix + "_tmp", 'w')
    try:
        p1 = Popen(cmd_args1, stdout=PIPE)
        p2 = Popen(cmd_args2, stdin=p1.stdout, stdout=tmp)
        p2.wait()
        tmp.close()
    finally:
        os.remove(bam)
        os.rename(bam_prefix + "_tmp", bam_prefix + "_frag_sort.bam")
        pysam.index(bam_prefix + "_frag_sort.bam")

def combineWigs(wig):
    out_wig = open(wig + ".wig", 'w')
    cmd_args = 'cat ' + wig + '/*'
    try:
        p = Popen(cmd_args, shell=True, stdout=out_wig)
        p.wait()
    finally:
        out_wig.close()
        shutil.rmtree(wig)
    
def calc_distance_score_worker(args):
    #pdb.set_trace()
    bam = args[0]
    sam_file = pysam.Samfile(bam, 'rb')
    resolution = args[1]
    curr_ref = args[2]
    curr_length = args[3]
    wig_prefix = args[4].split(".wig")[0]
    if not os.path.exists(wig_prefix): os.mkdir(wig_prefix)
    values = np.zeros(curr_length)
    distances = np.zeros(1000)
    distances_array_len = len(distances)
    res_count = 0
    #insert_size = 0
    read_len = 0
    sum_distances = 0
    num_reads = 0
    position = 0
    position_mean = 0
    position_var = 0
    sum_position_mean = 0
    sum_position_var = 0
    total_num_reads = 0
    value_ind = 0
    value_out = 0
    #first_pos = True
    #read_ind = 0
    wig_file = open("/".join([wig_prefix, curr_ref]), 'w')
    out = "fixedStep chrom={0} start=1 step={1} span={1}\n".format(curr_ref, resolution)
    wig_file.write(out)
    
    ## Create pileup iterator
    it = sam_file.pileup(reference=curr_ref)
    
    ## Loop through each position in iterator
    #print "Computing variances"
    for proxy in it:
        position = proxy.pos
        #if first_pos:
        value_ind = position / resolution
        #    first_pos = False
        #if res_count == 0: value_ind = value_ind + 1
        #assert(value_ind < curr_length)
        
        #pdb.set_trace()
        ## Loop through each read at current position
        for pread in proxy.pileups:
            #insert_size = pread.alignment.isize
            read_len = pread.alignment.alen
            #if insert_size > 50:
            if read_len > 50:
                #mid_distance = abs((pread.alignment.pos + (insert_size/2)) - position)
                mid_distance = abs((pread.alignment.pos + (read_len/2)) - position)
                distances[num_reads] = mid_distance
                #distances.append(mid_distance)
                #sum_distances = sum_distances + abs(read_mid)
                num_reads = num_reads + 1
                if num_reads == distances_array_len: break
            
        if num_reads > 10:
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
    
    max_value = np.max(values)
    scores = [0]
    
    if max_value > 0:
        scores = values
    #    scores = (max_value - values) / max_value
    
    #print "Writing scores"
    #for value in values:
    #    print>>wig_file, str(value)
    for score in scores:
        print>>wig_file, str(score)
        #wig_file.write(str(score) + "\n")
    wig_file.close()
    
def calc_distance_score(bam, wig, resolution, multi):
    #pdb.set_trace()
    sam_file = pysam.Samfile(bam, 'rb')
    #wig_file = open(wig, 'w')
    if not os.path.exists(wig): os.mkdir(wig)
    refs = sam_file.references
    ref_lengths = sam_file.lengths
    sam_file.close()
    
    pool = ""
    if multi:
        pool = multiprocessing.Pool(processes=1)
    #pdb.set_trace()
    for chr_index in range(len(refs)):
        #pdb.set_trace()
        curr_ref = refs[chr_index]
        curr_length = ref_lengths[chr_index]
        if resolution > 1:
            curr_length = curr_length / resolution
        arg = [bam, resolution, curr_ref, curr_length, wig]
        if multi:
            pool.apply_async(calc_distance_score_worker, (arg,))
        else:
            calc_distance_score_worker(arg)
            #pool.apply(calc_distance_score_worker, (arg,))
    if multi:
        pool.close()
        pool.join()
    
    combineWigs(wig)

    
def compute_mean_score(wig):
    N = 0
    sum_values = 0
    wig_file = open(wig)
    wig_file.next()
    for line in wig_file:
        value = float(line.split()[0])
        sum_values = sum_values + value
        N = N + 1
    wig_file.close()
    return sum_values / N

def compute_mean_score_set(set_name):
    path0 = "/".join([SIM_NUC_PATH, set_name])
    jitters = os.listdir(path0)
    jitters = [j for j in jitters if not re.search("desc", j)]
    #pdb.set_trace()
    jitters = [int(x) for x in jitters]
    
    jitters.sort()
    means = np.zeros((1000, len(jitters)))
    cind = 0
    for jitter in jitters:
        jitter = str(jitter)
        print jitter
        #pdb.set_trace()
        #if jitter == "desc": continue
        path1 = "/".join([path0, jitter])
        rind = 0
        files = os.listdir(path1)
        for file in files:
            if re.search("wig", file):
                means[rind, cind] = compute_mean_score("/".join([path1, file]))
                rind = rind + 1
        cind = cind + 1
    return(means)

def main(argv):
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', dest='bam')
    parser.add_argument('-w', dest='wig')
    parser.add_argument('-r', dest='resolution')
    args = parser.parse_args()
    
    
if __name__ == '__main__':
    main(sys.argv)
