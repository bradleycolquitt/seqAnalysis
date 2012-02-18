#!/usr/bin/env python

import sys, os
import pysam
import pdb

bed_dir = "/media/storage2/data/bed/"

def sam2bam(sam, bam):
    samfile = pysam.Samfile(sam, 'r')
    bamfile = pysam.Samfile(bam, 'wb', template=samfile)
    print "SAM -> BAM"
    for read in samfile:
        bamfile.write(read)
    bamfile.close()
    #os.remove(sam)
    
def proc(arg):
    bamfile = arg[0]
    rmdup = arg[1]
    #se = arg[2]
    print rmdup
    
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
            if not read.is_unmapped:
                mb.write(read)
        mb.close()
        mapped_fs = open(mapped_bam + "_stat", 'w')
        for line in pysam.flagstat(mapped_bam):
            mapped_fs.write(line)
        mapped_fs.close()    
    if not os.path.exists(rmdup_bam) and rmdup == "True":
        print "Removing duplicates..."
        pysam.rmdup("-S", mapped_bam, rmdup_bam)
        os.remove(mapped_bam)
        rmdup_fs = open(rmdup_bam, 'w')
        for line in pysam.flagstat(rmdup_bam):
            rmdup_fs.write(line)
        rmdup_fs.close()
        print "Sorting..."
        pysam.sort(rmdup_bam, sort_bam)
        os.remove(rmdup_bam)
    else:
        print "Sorting..."
        pysam.sort(mapped_bam, sort_bam)
        os.remove(mapped_bam)
    print "Indexing..."
    sort_bam = sort_bam + ".bam"
    pysam.index(sort_bam)
    #pysam.flagstat(bamfile)
    bamfile_fs = open(bamfile + "_stat", 'w')
    for line in pysam.flagstat(bamfile):
        bamfile_fs.write(line)
    bamfile_fs.close
    sort_bam_fs = open(sort_bam + "_stat", 'w')
    for line in pysam.flagstat(sort_bam):
        sort_bam_fs.write(line)
    sort_bam_fs.close()
    #os.remove(bamfile)

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
    
def main(argv):
    if argv[1] == "sam2bam":
        sam2bam(argv[2], argv[3])
    elif argv[1] == "proc":
        proc(argv[2:])
    elif argv[1] == "splitStrands":
        splitStrands(argv[2])
    elif argv[1] == "bam2bed":
        bam2bed(argv[2], )
    elif argv[1] == "removeWeirdChr":
        removeWeirdChr(argv[2])
        
if __name__ == "__main__":
    main(sys.argv)
