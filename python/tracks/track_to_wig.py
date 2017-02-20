#! /usr/bin/env python

import numpy as np
import genome.db
import genome.wig

import subprocess as sp

import pdb
import os
import shutil

def write_out_tracks_by_chrom(gdb, tmp_dir):
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)
    os.mkdir(tmp_dir)
    for chrom in gdb.get_chromosomes():
        vals = test_track.get_nparray(chrom.name)
        vals1 = vals.astype(np.float32)
        out = "/".join([tmp_dir, '%s.wig.gz' % chrom.name])
        genome.wig.write_float32(out, vals1, chrom.name)

def combine_chrom_wigs(tmp_dir, out_fname):
    out_file = open(out_fname, 'wb')
    out_files = os.listdir(tmp_dir)
    out_files = ["/".join([tmp_dir, f]) for f in out_files]
    args1 = ['zcat'] + out_files
    p1 = sp.Popen(args1, shell=False, stdout=sp.PIPE)
    args2 = ['gzip']
    p2 = sp.Popen(args2, shell=False, stdin=p1.stdout, stdout=out_file)
    p2.wait()

    shutil.rmtree(tmp_dir)
def convert_wig_to_bigwig(fname):
    out_fname = fname.strip('.wig.gz') + '.bw'
    args = ['wigToBigWig', fname, '/media/data2/assembly/lonStrDom1/chrom.sizes',]

def main():
    # 'connect' to the genome database
    gdb = genome.db.GenomeDB()

    # open data 'tracks' for DNase and MNase
    test_track = gdb.open_track('test.h5')
    tmp_dir = '/media/data2/hdf5/lonStrDom1/tmp'

    write_out_tracks_by_chrom(gdb, tmp_dir)

    print "Combining wigs to one file..."
    combine_chrom_wigs(tmp_dir, out_fname)
    convert_wig_to_bigwig(out_fname)

main()
