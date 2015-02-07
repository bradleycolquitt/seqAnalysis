#! /usr/bin/env python

#Read in filenames within directory
#Split fnames to get list of scaffolds
#For each scaffold, identify highest scoring chromosome
#return (scaffold, chrom, cstart, cend)


import sys
import os
import pandas as pd
from bx.align import lav
from operator import itemgetter
from fastahack import FastaHack

#Loop through all lav files containig scaffold name
#Find max score, return chrom, cstart, cend

def sum_scores(aligns):
    """ Sum scores and find minimum start and maximum end of lav file
    Args:
        aligns (lav.Reader record) alignment blocks(as structured by bx.align)

    Returns:
        tuple:
            summed scores of all alignment blocks
            minimum start position of query on target
            maximum end position of query on target
    """
    ss = 0
    min_start = 1E9
    max_end = 0

    for align in aligns:
        ss += align.score
        if align.components[0].start < min_start: min_start = align.components[0].start
        if align.components[0].get_end() > max_end: max_end = align.components[0].get_end()
    return (ss, min_start, max_end)

def find_highest_score(scaffold, fnames):
    """"""
    indices = [i for i, s in enumerate(fnames) if scaffold in s]
    fnames_wscaf = itemgetter(*indices)(fnames)

    sum_score = 0
    src = ""
    start = 0
    end = 0

    max_score = 0
    max_src = ""
    max_start = 0
    max_end = 0
    for fname in fnames_wscaf:
        lav_obj = lav.Reader(file(fname))
        try:
            aligns = [block for block in lav_obj]
            (sum_score, start, end) = sum_scores(aligns)

            if sum_score > max_score:
                max_score = sum_score
                max_src = lav_obj.seq1_src
                max_start = start
                max_end = end
        except AssertionError, e:
            continue

    return (max_src, max_start, max_end)

    # TODO: turn dict into 4 column dataframe
    # sort by chrom, start
    # create dict of chrom:fasta_file
    # foreach row
    #    read in lav
    #    extract query (seq2) sequence
    #    add 500 bp 'N'
    #    append to fasta file givenby chrom
def create_merged_seq(scaf_chrom, lav_dir):

    d = pd.DataFrame(scaf_chrom)
    return d


def main(argv):

    fnames = os.listdir(argv[1])
    fnames_scaf = [x.split("-")[1].split(".lav")[0] for x in fnames]
    fnames_scaf_set = set(fnames_scaf)

    scaf_chrom = {}
    for scaffold in fnames_scaf_set:
        print scaffold
        scaf_chrom[scaffold] = find_highest_score(scaffold, fnames)
    return scaf_chrom
    d = create_merged_seq(scaf_chrom)


if __name__ == "__main__":
    main(sys.argv)
