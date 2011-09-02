#!/usr/bin/env python

import sys, os
import pysam

bed_dir = "/media/storage2/data/bed/"

def sam2bam(sam, bam):
    samfile = pysam.Samfile(sam, 'r')
    bamfile = pysam.Samfile(bam, 'wb', template=samfile)
    print "SAM -> BAM"
    for read in samfile:
        bamfile.write(read)
    bamfile.close()
    #os.remove(sam)
    
def proc(bamfile):
    bam_dir = "/".join(bamfile.split("/")[:-1]) + "/"
    bam_prefix = os.path.basename(bamfile).split(".bam")[0]
    mapped_bam = bam_dir + bam_prefix + "_mapped.bam"
    rmdup_bam = bam_dir + bam_prefix + "_rmdup.bam"
    sort_bam = bam_dir + bam_prefix + "_sort"
    
    if not os.path.exists(mapped_bam):
        print "Removing unmapped..."
        bam = pysam.Samfile(bamfile, 'rb')
        mb = pysam.Samfile(mapped_bam, 'wb', template=bam)
        for read in bam:
            if read.is_unmapped or read.mate_is_unmapped:
                continue
            mb.write(read)
        mb.close()
    if not os.path.exists(rmdup_bam):
        
        print "Removing duplicates..."
        pysam.rmdup("-S", mapped_bam, rmdup_bam)
    #os.remove(mapped_bam)
    #print "Sorting..."
    pysam.sort(rmdup_bam, sort_bam)
    #os.remove(rmdup_bam)
    print "Indexing..."
    sort_bam = sort_bam + ".bam"
    pysam.index(sort_bam)
    #os.remove(bamfile)

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

def main(argv):
    if argv[1] == "sam2bam":
        sam2bam(argv[2], argv[3])
    elif argv[1] == "proc":
        proc(argv[2])
    elif argv[1] == "splitStrands":
        splitStrands(argv[2])
    elif argv[1] == "bam2bed":
        bam2bed(argv[2], )

if __name__ == "__main__":
    main(sys.argv)
