#!/usr/bin/env python

import os, sys
import argparse
import tables as tb
from subprocess import Popen

track_dir = '/media/storage2/data/h5'
lib_dir = '/home/user/lib/for_nuc'
pos_dir = '/media/storage2/data/nuc_pos'

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', dest="track_in", help="input track file")
    parser.add_argument('-b', dest="bed_set", help="bed set")
    parser.add_argument('-r', dest="reverse", help="number of bases to extend up from start")
    parser.add_argument('-f', dest="forward", help="number of bases to extend down from start")
                        
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
    
    out_path = "/".join([pos_dir, track_prefix, args.bed_set])
    if not os.path.exists(out_path): os.makedirs(out_path)
    
    for track in nodes:
        print track
        for bed_file in bed_files:
            out_file = "/".join([out_path, "_".join([track, bed_file, "R" + args.reverse, "F" + args.forward])]) 
            if os.path.exists(out_file): continue
            print "--", bed_file
            cmd_args = ['NucRegionPositions',
                '-i', "/".join([track_dir, args.track_in]),
                '-t', track,
                '-o', out_file,
                '-f', args.forward,
                '-r', args.reverse,
                '-b', "/".join([bed_path, bed_file])]
            try:
                p = Popen(cmd_args)
                p.wait()
            except:
                raise
    
if __name__ == '__main__':
    main(sys.argv)
