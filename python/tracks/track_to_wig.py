#! /usr/bin/env python

import pdb
import os
import shutil
import tempfile
import argparse

import numpy as np
import genome.db
import genome.wig

import subprocess as sp

def write_out_tracks_by_chrom(gdb, track, tmp_dir):
    for chrom in gdb.get_chromosomes():
        #if chrom.name != "scaffold_0":
        #    continue
        vals = track.get_nparray(chrom.name)
        vals1 = vals.astype(np.float32)
        out = "/".join([tmp_dir, '%s.wig.gz' % chrom.name])
        genome.wig.write_float32(out, vals1, chrom.name)

def combine_chrom_wigs(tmp_dir):
    #pdb.set_trace()
    out_files = os.listdir(tmp_dir)
    out_files = ["/".join([tmp_dir, f]) for f in out_files]

    out_fname = "/".join([tmp_dir, "combined.wig.gz"])
    out_file = open(out_fname, 'wb')
    args1 = ['zcat'] + out_files
    p1 = sp.Popen(args1, shell=False, stdout=sp.PIPE)
    args2 = ['gzip']
    p2 = sp.Popen(args2, shell=False, stdin=p1.stdout, stdout=out_file)
    p2.wait()
    return out_fname

    #shutil.rmtree(tmp_dir)
def convert_wig_to_bigwig(in_wig, out_bw, assembly):
    #out_fname = fname.strip('.wig.gz') + '.bw'
    args = ['wigToBigWig', in_wig, '/media/data2/assembly/%s/chrom.sizes' % assembly, out_bw]
    p = sp.Popen(args)
    p.wait()

def run(track_name, gdb, assembly):
    track = gdb.open_track(track_name)
    tmp_dir_prefix = '/media/data2/hdf5/lonStrDom1/tmp'
    tmp_dir = tempfile.mkdtemp(prefix=tmp_dir_prefix)

    bw_prefix = '/media/data2/wig'
    out_fname = "/".join([bw_prefix, assembly, track_name.strip(".h5") + ".bw"])
    write_out_tracks_by_chrom(gdb, track, tmp_dir)
    print "Combining wigs to one file..."
    combined_fname = combine_chrom_wigs(tmp_dir)
    print "Converting to bigwig..."
    convert_wig_to_bigwig(combined_fname, out_fname, assembly)
    shutil.rmtree(tmp_dir)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', dest='track')
    parser.add_argument('--assembly', help="assembly to use", default='lonStrDom1')
    args = parser.parse_args()




    gdb = genome.db.GenomeDB()
    run(args.track, gdb, args.assembly)
    # bw_prefix = '/media/data2/wig'
    # out_fname = "/".join([bw_prefix, args.assembly, args.track.strip(".h5") + ".bw"])

    # # open data 'tracks' for DNase and MNase
    # track = gdb.open_track(args.track)
    # tmp_dir_prefix = '/media/data2/hdf5/lonStrDom1/tmp'
    # tmp_dir = tempfile.mkdtemp(prefix=tmp_dir_prefix)

    # write_out_tracks_by_chrom(gdb, track, tmp_dir)
    # print "Combining wigs to one file..."
    # combined_fname = combine_chrom_wigs(tmp_dir)
    # print "Converting to bigwig..."
    # convert_wig_to_bigwig(combined_fname, out_fname, args.assembly)

    # shutil.rmtree(tmp_dir)

main()
