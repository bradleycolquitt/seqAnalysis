#!/usr/bin/env python

import sys, os, re
import pysam
import pdb
import numpy
import tables as tb
import track_util
from subprocess import Popen

bed_dir = "/media/storage2/data/bed/"

def sam2bam(sam, bam, errorlog):
    samfile = pysam.Samfile(sam, 'r')
    bamfile = pysam.Samfile(bam, 'wb', template=samfile)
    if errorlog == "stderr":
        errorlog = sys.stderr
    print>>errorlog, "SAM -> BAM"
    read_number = 0
    try:
        for read in samfile:
            read_number = read_number + 1
            bamfile.write(read)
        bamfile.close()
    except:
        raise("Error converting SAM to BAM: read {0}\n".format(read_number))
    
def proc(arg):
    bamfile = arg[0]
    rmdup = arg[1]
    errorlog = arg[2]
    if errorlog == "stderr":
        errorlog = sys.stderr
    if rmdup == "False": rmdup = False
    
    bam_dir = "/".join(bamfile.split("/")[:-1]) + "/"
    bam_prefix = os.path.basename(bamfile).split(".bam")[0]
    mapped_bam = bam_dir + bam_prefix + "_mapped.bam"
    rmdup_bam = bam_dir + bam_prefix + "_rmdup.bam"
    sort_bam = bam_dir + bam_prefix + "_sort"
    
    stat_dir = bam_dir + "stat/"
    if not os.path.exists(stat_dir): os.makedirs(stat_dir)
    
    if not os.path.exists(mapped_bam):
        print>>errorlog, "Removing unmapped..."
        mapped = 0
        unmapped = 0
        bam = pysam.Samfile(bamfile, 'rb')
        mb = pysam.Samfile(mapped_bam, 'wb', template=bam)
        try:
            for read in bam:
                if not read.is_unmapped:
                    mapped = mapped + 1
                    mb.write(read)
                else:
                    unmapped = unmapped + 1
        except:
            errorlog.write("Failed to remove unmapped reads: read number {0}\n".format(mapped + unmapped))
            raise
        else:
            errorlog.write("Unmapped read removal successful: Mapped {0}/Unmapped {1}\n".format(mapped, unmapped))

        bam.close()
        mb.close()
    
    if not os.path.exists(sort_bam + ".bam"):

        print>>errorlog, "Sorting..."

        try:
            cmd_args = ['java', '-Xmx2g', '-jar', '/seq/picard/SortSam.jar',
                        "=".join(["INPUT", mapped_bam]),
                        "=".join(["OUTPUT", sort_bam + ".bam"]),
                        "=".join(["SORT_ORDER", "coordinate"])]
            p = Popen(cmd_args, stdout=errorlog, stderr=errorlog)
            p.wait()
        except:
            errorlog.write("Sorting failed.\n")
            raise
        else:
            os.remove(mapped_bam)

    if not os.path.exists(rmdup_bam) and rmdup:   
        print "Removing duplicates..."
        rmdup_metrics = stat_dir + bam_prefix + "_rmdup_metrics"
        cmd_args = ['java', '-Xmx2g', '-jar', '/seq/picard/MarkDuplicates.jar',
                    "=".join(["INPUT", sort_bam + ".bam"]),
                    "=".join(["OUTPUT", rmdup_bam]),
                    "=".join(["METRICS_FILE", rmdup_metrics]),
                    "=".join(["REMOVE_DUPLICATES", "true"]),
                    "=".join(["ASSUME_SORTED", "true"])]
        try:
            p = Popen(cmd_args, stdout=errorlog, stderr=errorlog)
            p.wait()
        except:
            errorlog.write("Failed to remove duplicates.\n")
            raise
         
        try:
            print>>errorlog, "Indexing..."
            pysam.index(rmdup_bam)
        except SamtoolsError as detail:
            print>>errorlog, "Indexing failed: ",detail
    else:
        try:
            print>>errorlog, "Indexing..."
            sort_bam = sort_bam + ".bam"
            pysam.index(sort_bam)
        except SamtoolsError as detail:
            print>>errorlog, "Indexing failed: ", detail
    
 
    bamfile_fs = open(bam_dir + "stat/" + bam_prefix + "_stat", 'w')
    for line in pysam.flagstat(bamfile):
        bamfile_fs.write(line)
    bamfile_fs.close()

    return 0

# Takes a BAMtoWindow or similar object and given read
# Return midpoint as specified in obj
def read_mid_compute(obj, read):
    read_mid = 0
    
    if obj.ends:
        if read.is_reverse:
            read_mid = read.aend + 1
        else:
            read_mid = read.pos + 2
    elif obj.pe:
        if obj.extend > 0:
            if read.is_read1:
                if read.is_reverse:
                    read_mid = read.aend - (obj.extend / 2)
                else:
                    read_mid = read.pos + (obj.extend / 2)
        else:
            if read.is_paired and read.is_read1:
                if read.is_reverse:
                    read_mid = read.aend + read.isize / 2
                else:
                    read_mid = read.pos + read.isize / 2
            else: return -1
    else:
        if read.is_reverse:
            read_mid = read.aend - (obj.extend / 2)
        else:
            read_mid = read.pos + (obj.extend / 2)
    return read_mid

def filterISize(inbam, outbam, imin, imax):
    inbam_file = pysam.Samfile(inbam, "rb")
    outbam_file = pysam.Samfile(outbam, "wb", template=inbam_file)
    curr_ref = -1
    isize = 0
    for read in inbam_file.fetch():
        if read.tid != curr_ref:
            curr_ref = read.tid
            print curr_ref
        isize = abs(read.isize)   
        if isize > imin and isize < imax: 
            outbam_file.write(read)
    inbam_file.close()
    outbam_file.close()
    
def proc_sam(arg):
    samfile = arg[0]
    rmdup = arg[1]
    #se = arg[2]
    print samfile
    print rmdup
    
    sam_dir = "/".join(samfile.split("/")[:-1]) + "/"
    sam_prefix = os.path.basename(samfile).split(".sam")[0]
    mapped_sam = sam_dir + sam_prefix + "_mapped.sam"
    rmdup_sam = sam_dir + sam_prefix + "_rmdup.sam"
    sort_sam = sam_dir + sam_prefix + "_sort"
    
    if not os.path.exists(mapped_sam):
        print "Removing unmapped..."
        sam = pysam.Samfile(samfile, 'r')
        mb = pysam.Samfile(mapped_sam, 'w', template=sam)
        for read in sam:
            if not read.is_unmapped:
                mb.write(read)
        mb.close()
        print "Finished removing unmapped."
    if not os.path.exists(rmdup_sam) and rmdup == "True":
        print "Removing duplicates..."
        pysam.rmdup("-S", mapped_sam, rmdup_sam)
        os.remove(mapped_sam)
        print "Sorting..."
        pysam.sort(rmdup_sam, sort_sam)
        os.remove(rmdup_sam)
    else:
        print "Sorting..."
        pysam.sort(mapped_sam, sort_sam)
        os.remove(mapped_sam)
    print "Indexing..."
    sort_sam = sort_sam + ".sam"
    pysam.index(sort_sam)

    samfile_fs = open(samfile + "_stat", 'w')
    for line in pysam.flagstat(samfile):
        samfile_fs.write(line)
    samfile_fs.close
    sort_sam_fs = open(sort_sam + "_stat", 'w')
    for line in pysam.flagstat(sort_sam):
        sort_sam_fs.write(line)
    sort_sam_fs.close()


def splitStrands(bam):
    bam_base = bam.split(".bam")[0]
    plus_bam = bam_base + "_plus.bam"
    minus_bam = bam_base + "_minus.bam"
    samfile = pysam.Samfile(bam, 'rb')
    plus_file = pysam.Samfile(plus_bam, 'wb', template=samfile)
    minus_file = pysam.Samfile(minus_bam, 'wb', template=samfile)
    it = samfile.fetch()
    for read in it:
        if not read.is_reverse:
            plus_file.write(read)
        else:
            minus_file.write(read)
    samfile.close()
    plus_file.close()
    minus_file.close()
    
def bam2bed(sam, pe = True, region = None):
    samfile = pysam.Samfile(sam, 'rb')
    if region != None:
        it = samfile.fetch(region=region)
    else:
        it = samfile.fetch()
    # calculate the end position and print out BED
    take = (0, 2, 3) # CIGAR operation (M/match, D/del, N/ref_skip)
    bed = bed_dir + sam.split(".")[0] + ".bed"
    bedfile = open(bed, 'w')
    print "Writing BED"
    chrom = ""
    for read in it:
        tmp = samfile.getrname(read.rname)
        if tmp != chrom:
            print tmp
            chrom = tmp
        if read.is_unmapped or read.is_read2: continue
        # compute total length on reference
        t = 0
        if not pe:
            t = sum([ l for op,l in read.cigar if op in take ])
        else:
            if read.is_reverse:
                t = -1 * read.isize
                strand = "-"
            else :
                t = read.isize
                strand = "+"
        
        bedfile.write("%s\t%d\t%d\t%s\t%d\t%c\n" %\
                          ( samfile.getrname( read.rname ),
                            read.pos, read.pos+t, read.qname,
                            read.mapq, strand) )            
def removeWeirdChr(sam):
    ind = range(1,20)
    ind.append('X')
    ind.append('Y')
    chrs = ["".join(['chr', str(chr)]) for chr in ind]
    samfile = pysam.Samfile(sam, 'rb')
    head_in = samfile.header
    head_out = {}
    head_out['HD'] = head_in['HD']
    head_out['PG'] = head_in['PG']
    chrs_ind = []

    for index, ref in enumerate(head_in['SQ']):
        if ref['SN'] in chrs:
            chrs_ind.append(index)

    head_out['SQ'] = head_in['SQ'][0:22 ]

        
    out = os.path.basename(sam) + "_std"
    outfile = pysam.Samfile(out, 'wb', header=head_out)
    for chr in chrs:
        print chr
        it = samfile.fetch(reference=chr)
        for read in it:
            outfile.write(read)
            
    samfile.close()
    outfile.close()

def computeChrMean(bam):
    bamfile = pysam.Samfile(bam, 'rb')
    refs = bamfile.references
    ref_lengths = bamfile.lengths
    outfile = open(bam + ".chr_means_0", 'w')
    total = 0
    
    for index in range(len(refs)):
        ref_length = ref_lengths[index]
        ref = refs[index]
        print ref

        for i in range(int(3E6), int(3.1E6)):
            print i
            total += bamfile.count(ref, i, i+1)
            if (i % 1E7) == 0 : print "="
        pdb.set_trace()
        total /= (ref_length - 3E6)  
        outfile.write("\t".join([ref, str(total)]) + "\n")
                      

def extractInsertSizes(sam, wsize, output):
    sam_file = pysam.Samfile(sam, 'rb')
    wig_file = open(output, 'w')
    refs = sam_file.references
    ref_lengths = sam_file.lengths
    
    for chr_index in [0]:
        ## Extract and initialize variables
        curr_ref = refs[chr_index]
        print curr_ref
    
        curr_length = ref_lengths[chr_index]
        values = numpy.zeros(curr_length / wsize)
        sum_inserts = 0
        num_reads = 0
        value_ind = 0
        
        out = "fixedStep chrom={0} start=1 step={1} span={1}\n".format(curr_ref, wsize)
        wig_file.write(out)
    
        ## Create pileup iterator
        it = sam_file.pileup(reference=curr_ref)
        
        ## Loop through each position in iterator
        value_ind = 0
        value_ind_store = 0
        update = 1

        for proxy in it:
            assert(value_ind < curr_length)
            
            ## Loop through each read at current position
            for pread in proxy.pileups:
                sum_inserts = sum_inserts + abs(pread.alignment.isize)
                num_reads = num_reads + 1
                
            ## Compute average insert size at position
            value_ind = proxy.pos / wsize
            if value_ind_store == 0: value_ind_store = value_ind
            if value_ind == value_ind_store:
                update += 1
                values[value_ind] += float(sum_inserts) / float(num_reads)
                
            else:
                values[value_ind_store] /= float(update)
                if values[value_ind_store] > 1000: pdb.set_trace()
                update = 0
            value_ind_store = value_ind    
            
            sum_inserts = 0
            num_reads = 0
            
        
        values[value_ind_store] /= float(update)
        for value in values:
        
            wig_file.write(str(value) + "\n")
            
    wig_file.close()
    cmd_args = ['igvtools', 'tile', output, output + ".tdf", 'mm9']
    p = Popen(cmd_args)
    p.wait()

def main(argv):
    if argv[1] == "sam2bam":
        sam2bam(argv[2], argv[3], argv[4])
    elif argv[1] == "proc":
        proc(argv[2:])
    elif argv[1] == "splitStrands":
        splitStrands(argv[2])
    elif argv[1] == "bam2bed":
        bam2bed(argv[2], )
    elif argv[1] == "removeWeirdChr":
        removeWeirdChr(argv[2])
    elif argv[1] == "extractInsertSizes":
        extractInsertSizes(argv[2], argv[3], argv[4])
        
if __name__ == "__main__":
    main(sys.argv)
