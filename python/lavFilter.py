#! /usr/bin/env python

## First pass
# Identify target-query pair with highest match and write to tabbed file

## Second pass
# Identify longest alignment block for each pair and write out target-based start and stop for each query



import sys
import os
import pdb
import numpy as np
import pandas as pd
import multiprocessing as mp
from bx.align import lav
from operator import itemgetter
from fastahack import FastaHack

import matplotlib.pyplot as plt
from PIL import Image
from PIL import ImageDraw

#### Utilities
def plot_align(lav_file):
    lav_file.next()
    seq1_len = lav_file.seq1_end - lav_file.seq1_start
    seq2_len = lav_file.seq2_end - lav_file.seq2_start
    im = Image.new('RGBA', (seq1_len, seq2_len), (255,255,255,255))
    draw = ImageDraw.Draw(im)

    for align in lav_file:
        draw.line(((align.components[0].start, align.components[1].start), (align.components[0].get_end(), align.components[1].get_end())), fill=255, width=1)
    plt.imshow(np.asarray(im), origin='lower')
    plt.show()

def lav_to_df(lav_file, out):
    out_file = open(out, 'w')
    out_line = ""
    for align in lav_file:
        out_line = "\t".join(map(str, [align.components[0].start, align.components[0].get_end(), align.components[1].start, align.components[1].get_end()])) + "\n"
        out_file.write(out_line)

class LavProcess:
    def __init__(self, lav_dir):
        ### Read in filenames within lav directory
        self.lav_dir = os.path.abspath(lav_dir)
        self.fnames = [os.path.join(self.lav_dir, f) for f in  os.listdir(self.lav_dir)]

        ### Split fnames to get list of queries
        fnames_query = []
        for x in self.fnames:
            try:
                fnames_query.append(x.split("-")[1].split(".lav")[0])
            except:
                pass
        self.query_set = set(fnames_query)

        ### Output files
        self.out_dir = "/media/data2/assembly/lonStrDom1/lastz/lav_proccessing"
        if not os.path.exists(self.out_dir): os.mkdir(self.out_dir)

    def find_best_targets(self):
        """Find best query-target pair as determined by sum of alignment block scores
        Returns:
            best_target (dict): score, best target, query start on target, query stop on target
        Side effect:
            Writes best_target to tabbed text file
        """

        # PARALLEL
        pool = mp.Pool(processes=2)
        inputs = []
        for query in self.query_set:
            indices = [i for i, s in enumerate(self.fnames) if query in s]
            fnames_wquery = itemgetter(*indices)(self.fnames)

            # workaround for cases when only one fnames_wquery
            if isinstance(fnames_wquery, str): fnames_wquery = (fnames_wquery, )

            inputs.append((query, fnames_wquery))

        print "Finding query-target pairs..."
        best_target = dict(pool.map(_find_highest_score, inputs))

        # SERIES
        # best_target = {}
        # for query in self.query_set:
        #     print query
        #     indices = [i for i, s in enumerate(self.fnames) if query + ".lav" in s]
        #     fnames_wquery = itemgetter(*indices)(self.fnames)

        #     # workaround for cases when only one fnames_wquery
        #     if isinstance(fnames_wquery, str): fnames_wquery = (fnames_wquery, ) #

        #     best_target[query] = find_highest_score(query, fnames_wquery)

        ### Write out best query-target pairs
        print "Writing query-target pairs..."
        self.best_file = "/".join([self.out_dir, "lav_pairs.txt"])
        out = open(self.best_file, "w")
        for query,value in best_target.iteritems():
            out.write("\t".join([query] + map(str, list(value))) + "\n")
        out.close()

        return best_target

    def find_longest_alignments(self):
        """Groups alignment blocks into larget segments
        Returns:
            longest_alignments (dict):
                key - query name
                value (tuple) -
                    target-based start of longest segment
                    target-based end of longest segment
                    fraction coverage of longest segment over whole query
        """

        best_file = open(self.best_file)
        out = open(os.path.join(self.out_dir, "longest_segments.txt"), 'w')

        longest_alignments = {}
        for line in best_file:
            sline = line.split()
            fname = "{0}-{1}.lav".format(sline[2], sline[0])
            lav_file = lav.Reader(file(os.path.join(self.lav_dir, fname)))
            print fname

            a = lav_file.next()
            min_pos = a.components[0].start
            max_pos = a.components[0].get_end()
            query_length = lav_file.seq2_end - lav_file.seq2_start

            ### Group alignment blocks into larger segments
            last_target_pos = 0
            segments = []
            for align in lav_file:
                curr_target_pos = align.components[0].start
                if curr_target_pos > max_pos + 1000:
                    segments.append((min_pos, max_pos))
                    min_pos = curr_target_pos
                max_pos = align.components[0].get_end()

            ### Find longest segment
            lengths = [s[1] - s[0] for s in segments]
            longest = segments[lengths.index(max(lengths))]
            longest_alignments[sline[0]] = (longest[0], longest[1], float(longest[1] - longest[0]) / query_length)

            out.write("\t".join(sline + map(str, list(longest_alignments[sline[0]])) ) + "\n")
        out.close()
        return longest_alignments

def _find_highest_score(tup):
    query, fnames_query = tup
    return (query, find_highest_score(query, fnames_query))

def find_highest_score(query, fnames_wquery):
    sum_score = 0
    start = 0
    end = 0

    max_score = 0
    max_src = ""
    max_start = 0
    max_end = 0
    for fname in fnames_wquery:
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

    return (max_score, max_src, max_start, max_end)

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
    L = LavProcess(argv[1])

    #For each query, identify highest scoring target
    best_target = L.find_best_targets()
    #For each pair, find start/end coord on largest alignment block
    longest = L.find_longest_alignments()

    #d = create_merged_seq(scaf_chrom)


if __name__ == "__main__":
    main(sys.argv)
