#!/usr/bin/env python


import sys
import os
import re
import pdb
import argparse
import tables
import pysam

import numpy as np
import genome.db

MIN_MAP_QUAL = 10
MAX_VAL = 65535

def create_carray(track, chrom):
    atom = tables.UInt16Atom(dflt=0)

    zlib_filter = tables.Filters(complevel=1, complib="zlib")

    # create CArray for this chromosome
    shape = [chrom.length]
    carray = track.h5f.createCArray(track.h5f.root, chrom.name,
                                    atom, shape, filters=zlib_filter)

    return carray




def get_sam_iter(samfile, chrom):
    try:
        sam_iter = samfile.fetch(reference=chrom.name,
                                 start=1, end=chrom.length)
    except ValueError:
        # could not find chromosome, try stripping leading 'chr'
        # E.g. for drosophila, sometimes 'chr2L' is used but
        # othertimes just '2L' is used. Annoying!
        chrom_name = chrom.name.replace("chr", "")
        sys.stderr.write("WARNING: %s does not exist in BAM file, "
                         "trying %s instead\n" % (chrom.name, chrom_name))

        try:
            sam_iter = samfile.fetch(reference=chrom_name,
                                     start=1, end=chrom.length)
        except ValueError:
            sys.stderr.write("WARNING: %s does not exist in BAM file, "
                             "skipping chromosome\n" % chrom_name)
            sam_iter = iter([])

    return sam_iter


def add_read_counts(chrom, fwd_array, rev_array, bam_filename):
    pdb.set_trace()
    samfile = pysam.Samfile(bam_filename, "rb")

    count = 0

    for read in get_sam_iter(samfile, chrom):
        count += 1
        if count > 100000:
            sys.stderr.write(".")
            count = 0

        if read.is_unmapped:
            continue

        if read.mapq < MIN_MAP_QUAL:
            # read has poor mapping quality
            continue

        # remember pysam pos starts at 0, not 1
        start = read.pos + 1
        end = start + read.alen - 1

        if end > chrom.length:
            sys.stderr.write("  warning: read coordinates %d-%d fall outside "
                             "%s range 1-%d" % (start, end, chrom.name,
                                                 chrom.length))
            end = chrom.length

        if start < 1:
            sys.stderr.write("  warning: read coordinates %d-%d fall outside "
                             "%s range 1-%d" % (start, end, chrom.name,
                                                 chrom.length))
            start = 1

        if read.is_reverse:
            rev_array[end-1] += 1
        else:
            fwd_array[start-1] += 1


def run_script(bam_filename, assembly):
    gdb = genome.db.GenomeDB(assembly=assembly)

    bam_base = os.path.basename(bam_filename).strip(".bam")
    fwd_track_name = bam_base + "-fwd"
    rev_track_name = bam_base + "-rev"

    gdb.delete_track(fwd_track_name)
    gdb.delete_track(rev_track_name)

    fwd_track = gdb.create_track(fwd_track_name)
    rev_track = gdb.create_track(rev_track_name)

    for chrom in gdb.get_all_chromosomes():
        sys.stderr.write("%s\n" % chrom.name)

        fwd_carray = create_carray(fwd_track, chrom)
        rev_carray = create_carray(rev_track, chrom)

        # use uint32 even though we will store values as uint16s
        # to allow for possible overflows
        fwd_array = np.zeros(chrom.length, np.uint32)
        rev_array = np.zeros(chrom.length, np.uint32)

        # fill with values

        sys.stderr.write("  %s\n  " % bam_filename)
        add_read_counts(chrom, fwd_array, rev_array, bam_filename)
        sys.stderr.write("\n")

        # threshold values to avoid integer overflow when we store them
        threshold_large_vals(fwd_array)
        if rev_track:
            threshold_large_vals(rev_array)

        fwd_carray[:] = fwd_array
        rev_carray[:] = rev_array

    fwd_track.close()
    rev_track.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', dest='in_dir')

    parser.add_argument('--assembly', help="assembly to use", default=None)
    args = parser.parse_args()

    full_in_dir = "/".join([os.getcwd(), args.in_dir])
    in_files = os.listdir(full_in_dir)
    in_files = [f for f in in_files if re.search("bam$", f)]
    in_files = ["/".join([full_in_dir, f]) for f in in_files]

    for f in in_files:
        run_script(f, args.assembly)


if __name__ == '__main__':
    main()