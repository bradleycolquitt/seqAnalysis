#!/usr/bin/env python

import os, sys
import argparse
import tables as tb
from subprocess import Popen

track_dir = '/media/storage2/data/h5'
lib_dir = '/home/user/lib/for_nuc'
phase_dir = '/media/storage2/analysis/nuc/phasogram'

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', dest="track_in", help="input track file")
    parser.add_argument('-b', dest="bed_set", help="bed set")
    parser.add_argument('-n', dest="thresh", help="pileup threshold")
    parser.add_argument('-l', dest="length", help="histogram length")
    
    args = parser.parse_args()
    
    track_path = "/".join([track_dir, args.track_in])
    track_prefix = args.track_in.split("_")[0]
    trk = tb.openFile(track_path)
    nodes = []
    for node in trk.listNodes("/"):
        nodes.append(node._v_name)
    trk.close()
    
    bed_path = "/".join([lib_dir, args.bed_set])
    bed_files = os.listdir(bed_path)
    
    out_path = "/".join([phase_dir, args.bed_set])
    if not os.path.exists(out_path): os.makedirs(out_path)
    
    for track in nodes:
        print track
        for bed_file in bed_files:
            out_file = "/".join([out_path, "_".join([track, bed_file, "N" + args.thresh, "L" + args.length])]) 
            if os.path.exists(out_file): continue
            print "--", bed_file
            cmd_args = ['Phasogram',
                '-i', "/".join([track_dir, args.track_in]),
                '-t', track,
                '-o', out_file,
                '-n', args.thresh,
                '-l', args.length,
                '-b', "/".join([bed_path, bed_file])]
            try:
                p = Popen(cmd_args)
                p.wait()
            except:
                raise
    
if __name__ == '__main__':
    main(sys.argv)
