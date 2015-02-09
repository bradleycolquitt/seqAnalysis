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
import bx_lav as lav
import fasta_utils as futil
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

#### Main class, functions
class LavProcess:
    def __init__(self, lav_dir):
        ### Read in filenames within lav directory
        self.lav_dir = os.path.abspath(lav_dir)
        self.fnames = [os.path.join(self.lav_dir, f) for f in  os.listdir(self.lav_dir)]
        self.query_fas = {}

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

        self.best_file = "/".join([self.out_dir, "lav_pairs.txt"])

        dec = "y"
        if os.path.exists(self.best_file):
            dec = raw_input("lav_pairs exists. Overwrite? [y/n]")
            if dec == "n": return

        # PARALLEL
        pool = mp.Pool(processes=10)
        inputs = []
        for query in self.query_set:
            indices = [i for i, s in enumerate(self.fnames) if query in s]
            fnames_wquery = itemgetter(*indices)(self.fnames)

            # workaround for cases when only one fnames_wquery
            if isinstance(fnames_wquery, str): fnames_wquery = (fnames_wquery, )

            inputs.append((query, fnames_wquery))

        print "Finding query-target pairs..."
        best_target = dict(pool.map(_find_highest_score, inputs))

        # # SERIES
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
            if value[0] == 0: continue
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
        skipped = open(os.path.join(self.out_dir, "skipped_pairs.txt"), 'w')
        out.write("\t".join(["query_name",
                             "score",
                             "target_name",
                             "target_initial_start",
                             "target_initial_end",
                             "target_final_start",
                             "target_final_end",
                             "query_strand",
                             "query_coverage"]) + "\n")
        longest_alignments = {}
        for line in best_file:
            sline = line.split()
            fname = "{0}-{1}.lav".format(sline[2], sline[0])
            print fname
            score = int(sline[1])

            segments = []
            #span = []

            ### Initial
            lav_file = lav.Reader(file(os.path.join(self.lav_dir, fname)))
            a = lav_file.next()
            self.query_fas[lav_file.seq2_src] = lav_file.seq2_filename

            min_pos = a.components[0].start
            max_pos = a.components[0].get_end()
            query_length = lav_file.seq2_end - lav_file.seq2_start

            ### If poor coverage, skip
            if (score / 100 < query_length * 0.2):
                print "Low score: {0}".format(fname)
                skipped.write("\t".join(sline) + "\tlow_score\n")
                continue

            ### Initial segment read and merge
            gap = 1000
            segments.append([])
            strand = "+"

            for align in lav_file:
                curr_target_pos = align.components[0].start
                if abs(curr_target_pos - max_pos) >  gap:
                    segments[0].append((min_pos, max_pos, align.components[1].strand))
                    min_pos = curr_target_pos
                max_pos = align.components[0].get_end()
                strand = align.components[1].strand

            segments[0].append((min_pos, max_pos, strand))
            segments[0] = sorted(segments[0], key=lambda x: x[0])

            ### Drop short segments
            segments[0] = [s for s in segments[0] if float(s[1]-s[0]) / query_length > 0.001]
            #span.append(segments[0][-1][1] - segments[0][0][0])
            #pdb.set_trace()
            try:
                final = self.merge_segments(segments, query_length)
            except:
                skipped.write("\t".join(sline) + "\tbad_match\n")
                continue

            if len(final) > 1:
                spans = [x[1] - x[0] for x in final]
                final = final[spans.index(max(spans))]
            else:
                final = final[0]

            coverage = round(float(final[1] - final[0]) / query_length, 2)
            longest_alignments[sline[0]] = (final[0], \
                                            final[1], \
                                            final[2], \
                                            coverage)
            out.write("\t".join(sline + map(str, list(longest_alignments[sline[0]])) ) + "\n")

        out.close()
        return longest_alignments

    def merge_segments(self, segments, query_length):
        gap_range = [20000, 50000, 100000, 1000000]
        for i in range(1, len(gap_range) + 1):
            segment_iter = iter(segments[i-1])
            segment = next(segment_iter)

            min_pos = segment[0]
            max_pos = segment[1]
            strand = segment[2]
            segments.append([])

            for segment in segment_iter:
                curr_target_pos = segment[0]
                if curr_target_pos > max_pos + gap_range[i-1]:
                    segments[i].append((min_pos, max_pos, strand))
                    min_pos = curr_target_pos
                    strand = segment[2]
                max_pos = segment[1]
            segments[i].append((min_pos, max_pos, strand))

            ### Drop short segments
            segments[i] = [s for s in segments[i] if float(s[1]-s[0]) / query_length > 0.01]

        return segments[-1]

    def create_merged_seq(self):
        """Writes merged target-ordered query sequence"""
        align = pd.read_table(os.path.join(self.out_dir, "longest_segments.txt"))
        align.sort(columns=["target_name", "target_final_start"], inplace=True)

        out = open(os.path.join(self.out_dir, "merged_seq.fa"), "w")

        curr_target_name = ""
        i = 0
        for row in align.iterrows():
            query = row[1]['query_name']
            target = row[1]['target_name']
            fa = FastaHack(self.query_fas[query])

            if curr_target_name != target:
                if i > 0:
                    out.write("\n")
                    i = 1
                out.write(">{0}\n".format(target))
                curr_target_name = target

            fasub = fa.get_sequence(query)
            if row[1]['query_strand'] == "-":
                fasub = futil.reverse_complement(fasub)
            fasub1 = futil.format_fasta(fasub)
            [out.write(x) for x in fasub1]

            ns = "N" * 10000
            ns1 = futil.format_fasta(ns)
            [out.write(x) for x in ns1]

        out.write("\n")

def _find_highest_score(tup):
    query, fnames_query = tup
    print query
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

def main(argv):
    L = LavProcess(argv[1])

    #For each query, identify highest scoring target
    best_target = L.find_best_targets()
    #For each pair, find start/end coord on largest alignment block
    longest = L.find_longest_alignments()

    print "Writing sequence"
    seq = L.create_merged_seq()


if __name__ == "__main__":
    main(sys.argv)
