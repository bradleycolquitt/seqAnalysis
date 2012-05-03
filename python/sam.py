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
    #pdb.set_trace()
    ## Remove poorly formatted records
    #errorlog.write("Removing poorly formatted reads...\n")
    #sam_in = open(sam, 'r')
    #log_path = "/".join([os.path.dirname(sam), "log", "sam2bam"])
    #if not os.path.exists(log_path): os.makedirs(log_path)
    #sam_error = open("/".join([log_path, os.path.basename(sam)]), 'a')
    #sam_out_name = sam + "_tmp"
    #sam_out = open(sam_out_name, 'w')
    #header_flag = re.compile("@")
    #sline = []
    #cigar_len = 0
    #seq_len = 0
    #line_number = 0
    #try:
    #    for line in sam_in:
    #        line_number = line_number + 1
    #        if not header_flag.search(line):
    #            sline = line.split()
    #            #if len(sline) < 14:
    #            #    sam_error.write(line)
    #            #    continue
    #            if sline[5] != '*':
    #                cigar_len = int(sline[5].split("M")[0])
    #                seq_len = len(sline[9])
    #                if cigar_len != seq_len: 
    #                    #pdb.set_trace()
    #                    sam_error.write(line)
    #                    continue
    #           
    #                
    #        sam_out.write(line)
    #except:
    #    error_log.write(">> Error cleaning SAM: line {0}\n".format(line_number))
    #else:
    #    error_log.write(">> Successfully cleaned SAM.\n")
    #finally:
    #    sam_in.close()
    #    sam_error.close()
    #    sam_out.close()
    ##os.rename(sam_out_name, sam)
    
    #samfile = pysam.Samfile(sam_out_name, "r")
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
    else:
        os.remove(sam)
        #os.remove(sam_out_name)
    
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
        #pdb.set_trace()
        print>>errorlog, "Sorting..."
        #cmd_args = ['samtools', 'sort', mapped_bam, sort_bam]
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
        #pysam.sort(rmdup_bam, sort_bam)
    if not os.path.exists(rmdup_bam) and rmdup:   
        print "Removing duplicates..."
        #cmd_args = ['samtools', 'rmdup', mapped_bam, rmdup_bam]
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
    
    #sort_bam_fs = open(sort_bam + "_stat", 'w')
    #for line in pysam.flagstat(sort_bam):
    #    sort_bam_fs.write(line)
    #sort_bam_fs.close()
    #os.remove(bamfile)
    return 0

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
    #pysam.flagstat(samfile)
    samfile_fs = open(samfile + "_stat", 'w')
    for line in pysam.flagstat(samfile):
        samfile_fs.write(line)
    samfile_fs.close
    sort_sam_fs = open(sort_sam + "_stat", 'w')
    for line in pysam.flagstat(sort_sam):
        sort_sam_fs.write(line)
    sort_sam_fs.close()
    #os.remove(samfile)

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
    #pdb.set_trace()
    for index, ref in enumerate(head_in['SQ']):
        if ref['SN'] in chrs:
            chrs_ind.append(index)
    #pdb.set_trace() 
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

def extractInsertSizes(sam, output, track_name):
    sam_file = pysam.Samfile(sam, 'rb')
    #h5_file = tb.openFile(output, 'a')
    wig_file = open(output, 'w')
    refs = sam_file.references
    ref_lengths = sam_file.lengths
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
        values = numpy.zeros(curr_length)
        sum_inserts = 0
        num_reads = 0
        value_ind = 0
        
        out = "fixedStep chrom={0} start=1 step={1} span={1}\n".format(curr_ref, 1)
        wig_file.write(out)
        ## Create pileup iterator
        it = sam_file.pileup(reference=curr_ref)
        
        ## Loop through each position in iterator
        #print "Computing"
        for proxy in it:
            assert(value_ind < curr_length)
            
            ## Loop through each read at current position
            for pread in proxy.pileups:
                sum_inserts = sum_inserts + abs(pread.alignment.isize)
                num_reads = num_reads + 1
                
            ## Computre average insert size at position
            #values[value_ind] = float(sum_inserts) / float(num_reads)
            value = float(sum_inserts) / float(num_reads)
            wig_file.write(str(value) + "\n")
            sum_inserts = 0
            num_reads = 0
            value_ind = value_ind + 1
            
        ## Write HDF5 track
        #print "Writing"
        #wig_file.write()
        #track_util.writeTrack(h5_file, track_name, curr_ref, values, 1)
    cmd_args = ['igvtools', 'tile', output, output + ".tdf", 'mm9']
    p = Popen(cmd_args)
    p.wait()
    #h5_file.flush()
    #h5_file.close()

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
