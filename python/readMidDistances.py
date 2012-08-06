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
from multiprocessing import Pool
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
    def write(self, total_reads, B, jitter_range, insert_size, jitter_advance, comment):
        self.file.write("Simulated nucleosomes\n")
        self.file.write("Total number of reads: " + str(total_reads) + "\n")
        self.file.write("Number of simulations: " + str(B) + "\n")
        self.file.write("Range of midpoint jitterle factors: " + str(" ".join(str(jitter_range))) + "\n")
        self.file.write("Standard deviation of insert sizes (mean 145): " + str(insert_size) + "\n")
        if jitter_advance > 0:
            self.file.write("Nucleosome advanced by draw from normal distribution with mean 300, sd %s\n" % str(jitter_advance))
        if comment != "none":
            self.file.write(comment + "\n")
            
# Simulate B multiple nucleosome BAMs over given set of jitters (jitter_range)         
def simulate_nuc_range_multiple(set_no, jitter_range=(8,16,32,64), total_reads=10000,
                                 B=100, insert_sd=10, advance=300, jitter_advance=0, comment="none"):

    read_path = "/".join([SIM_NUC_PATH, "set" + str(set_no)])
    if not os.path.exists(read_path): os.makedirs(read_path)
    else:
        dec = raw_input("set exists. Overwrite [y/n]? ")
        if dec == "n": return
    desc = desc_file(read_path)
    desc.write(total_reads, B, jitter_range, insert_sd, jitter_advance, comment)
    
    for jitter in jitter_range:
        print jitter
        jitter_path = "/".join([read_path, str(jitter)])
        if not os.path.exists(jitter_path): os.mkdir(jitter_path)
        pool = multiprocessing.Pool(processes=4)    
        args = []
        for b in xrange(B):
            out_bam = "/".join([jitter_path, "sim" + str(b)])
            args = [total_reads, jitter, insert_sd, advance, jitter_advance, out_bam]
            #pool.apply_async(sim_measure_worker, (args,))
            #pool.apply(sim_measure_worker, (args,))
            sim_measure_worker(args)
        pool.close()
        pool.join()
        
def sim_measure_worker(args):
    sim = simulate_nuc(args[0], args[1], args[2], args[3], args[4], args[5])
    #print "Sim run"
    sim.run()
    #print "Extend"
    #pdb.set_trace()
    extend_bam(args[5] + ".bam")
    #print "Distance"
    calc_distance_score(args[5] + "_ex.bam", args[5] + "_ex", 1, False)
          
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

# Take bed of queryname sorted paired end reads
# Output bed records with start = start of read1 , end = end of read2
# Workaround of failing -bedpe flag of bamToBed 
def stitchBed(infile, outfile):
    #pdb.set_trace()
    first = ""
    for line in infile:
        line = line.split()
        id = line[3].split("/")
        if id[1] == "1": first = {id[0]: line}
        elif id[1] == "2":
            if first.keys()[0] != id[0]: continue
            ind = (1, 2)
            if line[5] == "+": ind = (2, 1)
            out = "\t".join([line[0], first[id[0]][ind[0]], line[ind[1]], id[0], line[4], "+"]) + "\n"
            outfile.write(out)

def bamToFragmentBed(infile, outfile, size):
    out = ""
    start = 0
    end = 0
    if size == 0:
        for read in infile:    
            if read.is_read1:
                if read.is_reverse and read.isize < 0:
                    start = read.pos + read.isize
                    if start < 0: continue
                    end = read.pos
                    strand = "-"
                elif not read.is_reverse and read.isize > 0:
                    end = read.pos + read.isize
                    start = read.pos
                    strand = "+"
                else: continue
                if start >= end: continue
                out = "\t".join([infile.getrname(read.tid), str(start), str(end),
                                 read.qname, '0', strand]) + "\n"
                outfile.write(out)
    else:
        for read in infile:
            if read.tid < 0: continue
            if read.is_reverse:
                start = read.pos - size
                if start < 0: continue
                end = read.pos
                strand = "-"
            elif not read.is_reverse:
                end = read.pos + size
                start = read.pos
                strand = "+"
            else: continue
            if start >= end: continue
            out = "\t".join([infile.getrname(read.tid), str(start), str(end),
                            read.qname, '0', strand]) + "\n"
            outfile.write(out)

def bedToBam(bed_name, bam_name):
    try:
        out_bam = open(bam_name, 'w')
        print " ".join(["Converting", bed_name, "to", bam_name]) + "..."
        cmd_args = ['bedToBam', '-i', bed_name, '-g', '/seq/lib/mouse.mm9.genome_norand']
        p1 = Popen(cmd_args, stdout=out_bam)
        p1.wait()
    except:
        print "Failed to convert BED to BAM."
        out_bam.close()
        return
    else:
        print "BED to BAM successful."
        os.remove(bed_name)
        
# Convert paired bed BAM to BAM with extended reads
# In effect, convert pairs of reads to fragments of given insert size
# Use to prep for calc_distance_score

def trimToDyad_worker(args):
    bam = args[0]
    ref = args[1]
    it = bam.fetch(reference=ref)
    bed = args[2]
    surround = args[3]
    
    #aend is aligned end of read, 3' most end for both directions
    #pos is 0-based leftmost coordinate
    for read in it:
        if read.is_read1:
            if read.is_reverse and read.isize < 0:
                #pdb.set_trace()
                mid = read.aend + (read.isize / 2) 
                if mid < 0: continue
                strand = "-"
            elif not read.is_reverse and read.isize > 0:
                #pdb.set_trace()
                mid = read.pos + (read.isize / 2)  # 3 is correction for some weird behavior in bedToBam
                strand = "+"
            else: continue
            start = mid - surround
            end = mid + surround
            out = "\t".join([bam.getrname(read.tid), str(start), str(end),
                             read.qname, '0', strand]) + "\n"
            bed.write(out)
            
def trimToDyad(bam, bed, surround):
    out = ""
    mid = 0
    start = 0
    end = 0
    strand = ""
#    pdb.set_trace()

    refs = bam.references
    for ref in refs:
        args = [bam, ref, bed, surround]

        #pool.apply_async(trimToDyad_worker, ([it, bed, surround],))
        #pool.apply(trimToDyad_worker, (args,))
        trimToDyad_worker([bam, ref, bed, surround])
  
def extend_bam(bam, type, reheader, size=0):
    
    bam_prefix = bam.split(".bam")[0]
    bam_file = pysam.Samfile(bam, 'rb')
    tmp_name = bam_prefix + ".bed"
    tmp_bed = open(tmp_name, 'w')
    out_name = "_".join([bam_prefix, type]) + ".bam"
    out_bam = open(out_name, 'w')
    #pdb.set_trace()
    ## Convert BAM to temporary BED
    try:
        print "BAM to BED..."
        if type=="extend":
            bamToFragmentBed(bam_file, tmp_bed, size)
        elif type=="dyad":
            trimToDyad(bam_file, tmp_bed, 20)
    except:
        print "BAM to BED conversion failed."
        print ">> " + ":".join(sys.exc_info()[1])
        tmp_bed.close()
        out_bam.close()
        os.remove(tmp_name)
        return
    else:
        print "BAM to BED conversion successful."
        tmp_bed.close()
        #out_bam.close()
    
    ## Convert tmp bed to bam
    bedToBam(tmp_name, out_name)
    
    ## Replace header
    if reheader:
        cmd_args1 = ['samtools', 'view', '-h', bam]
        cmd_args2 = ['samtools', 'reheader', '-', out_name]
        tmp_name = bam_prefix + "_tmp"
        tmp = open(tmp_name, 'w')
        try:
            print "Reheader..."
            p1 = Popen(cmd_args1, stdout=PIPE)
            p2 = Popen(cmd_args2, stdin=p1.stdout, stdout=tmp)
            p2.wait()
        except:
            print "Failed reheader"
            tmp.close()
            os.remove(tmp_name)
            return
        else:
            
            #os.remove(bam)
            tmp.close()
            #os.rename(tmp_name, out_name)
    print "Sorting..."    
    pysam.sort(out_name, out_name + "_sort")
    os.rename(out_name + "_sort.bam", out_name)
    pysam.index(out_name)


        
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
    curr_length = args[3] / resolution
    wig_prefix = args[4].split(".wig")[0]
    report = args[5]
    
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
    position_means = 0
    position_vars = 0
    if resolution > 1:
        position_vars = np.zeros(resolution)
    total_num_reads = 0
    region_read_count = 0
    value_ind = 0
    value_out = 0
    #first_pos = True
    #read_ind = 0
    wig_file = 0
    if not report:
        wig_file = open("/".join([wig_prefix, curr_ref]), 'w')
        out = "fixedStep chrom={0} start=1 step={1} span={1}\n".format(curr_ref, resolution)
        wig_file.write(out)
    #pdb.set_trace()
    ## Create pileup iterator
    it = 0
    if report:
        it = sam_file.pileup(region=curr_ref)
        region_read_count = sam_file.count(region=curr_ref)
    else:
        it = sam_file.pileup(reference=curr_ref)
        region_read_count = sam_file.count(reference=curr_ref)
    positions = 0
    ## Loop through each position in iterator
    for proxy in it:
       # pdb.set_trace()
        #positions = positions + 1
        position = proxy.pos
        if not report:
            value_ind = position / resolution
       
        ## Loop through each read at current position
        for pread in proxy.pileups:
            read_len = pread.alignment.alen
            if read_len > 20:
                mid_distance = abs((pread.alignment.pos + (read_len/2)) - position)
                distances[num_reads] = mid_distance
                num_reads = num_reads + 1
                if num_reads == distances_array_len: break
            
        if num_reads > 10:
            #pdb.set_trace()
            position_mean = np.mean(distances[:num_reads])
            #position_var = sum(pow((distances[:num_reads] - position_mean), 2)) / (num_reads - 1)
            #position_var = np.sqrt(sum(pow((distances[:num_reads] - position_mean), 2)) / (num_reads - 1))
            #position_var = np.sqrt(np.var(distances[:num_reads]))
            position_var = np.std(distances[:num_reads], ddof=1)
            #total_num_reads = total_num_reads + num_reads
        #sum_position_mean = sum_position_mean + position_mean
        #sum_position_var = sum_position_var + position_var
        if resolution > 1:
            position_means[res_count] = position_mean
            position_vars[res_count] = position_var
        else:
            position_means = position_mean
            position_vars = position_var
        
        num_reads = 0
        position_mean = 0
        position_var = 0
        res_count = res_count + 1   
        if res_count == resolution:
            if resolution > 1:
                sum_position_mean = np.sum(position_means)
                sum_position_var = np.sum(position_vars)
            else:
                sum_position_mean = position_means
                sum_position_var = position_vars
                
            if sum_position_var > 0 and sum_position_mean <= 100:
                #value_out = (sum_position_mean / sum_position_var) / resolution
                #value_out = 1 / ((1 + sum_position_mean) * sum_position_var * resolution) 
                value_out = (100 - abs(sum_position_mean)) / (sum_position_var * resolution)
                #value_out = total_num_reads / (sum_position_var * resolution)
                #value_out = sum_position_mean / resolution
                #value_out = sum_position_var / resolution
            values[value_ind] = value_out
            if report:
                value_ind = value_ind + 1
            value_out = 0
            res_count = 0
            sum_position_mean = 0
            sum_position_var = 0
            total_num_reads = 0       
        #value_ind = value_ind + 1
        if value_ind == len(values):
            break
    scores = [0]
    max_value = np.max(values)
    #print positions
    if max_value > 0:
        #pdb.set_trace()
        #scores = values / region_read_count
        scores = values
    #pdb.set_trace()
    if report:
        return scores
    else:
        for score in scores:
            print>>wig_file, str(score)
        wig_file.close()
    
def calc_distance_score(bam, wig, resolution, multi):
    #pdb.set_trace()
    sam_file = pysam.Samfile(bam, 'rb')
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
        arg = [bam, resolution, curr_ref, curr_length, wig, False]
        if multi:
            pool.apply_async(calc_distance_score_worker, (arg,))
        else:
            calc_distance_score_worker(arg)
            #pool.apply(calc_distance_score_worker, (arg,)
        
    if multi:
        pool.close()
        pool.join()
    
    combineWigs(wig)

def calc_distance_score_feature(feature, bam, resolution, multi):
    feature_path = "/".join(["/home/user/lib/features_general", feature])
    feature_file = open(feature_path)
    sam_file = pysam.Samfile(bam, 'rb')
    out_path = "/".join(["/home/user/data/nuc/cv2", os.path.basename(bam)])
    if not os.path.exists(out_path): os.makedirs(out_path)
    out_file = open("/".join([out_path, feature]), 'w')
    
    pool = ""
    if multi:
        pool = multiprocessing.Pool(processes=4)
    #num_features = 0
    
    for feature in feature_file:
        #num_features = num_features + 1
        #pdb.set_trace()
        feature_split = feature.split()
        curr_length = int(feature_split[2]) - int(feature_split[1])
        curr_region = "{0}:{1}-{2}".format(feature_split[0], feature_split[1], feature_split[2])
        if resolution > 1:
            curr_length = curr_length / resolution
        arg = [bam, resolution, curr_region, curr_length, "none", True, True]
        scores = 0
        #pdb.set_trace()
        if multi:
            try: 
                scores = pool.apply_async(calc_distance_score_worker, (arg,)).get()
            except IndexError as e:
                print "Index Error: " + " ".join(feature_split)
            except ValueError as e:
                print "Value Error: " + " ".join(feature_split)
            #scores = pool.apply(calc_distance_score_worker, (arg,))
        else:
            scores = calc_distance_score_worker(arg)
            
        score_mean = np.mean(scores)
        out = "\t".join([feature.strip(), str(score_mean)]) + "\n"
        out_file.write(out)
   
        
    if multi:
        pool.close()
        pool.join()

def compute_mean_score(wig):
    N = 0
    sum_values = 0
    wig_file = open(wig)
    wig_file.next()
    for line in wig_file:
        value = float(line.split()[0])
        if value == 0: continue
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
    #pbd.set_trace()
    n_sim = len([f for f in os.listdir("/".join([path0, str(jitters[0])])) if re.search("wig", f)])
    means = np.zeros((n_sim, len(jitters)), dtype=np.float32) 
    cind = 0
    for jitter in jitters:
        jitter = str(jitter)
        print jitter
        #pdb.set_trace()
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
