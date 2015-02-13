#!/usr/bin/env python

import sys
import argparse
import tables as tb
import numpy as np
import re
import pdb
import track_util as tutil
import signal_utils

class trackProcess:
    def __init__(self, in_track, out_track, smooth):
        self.in_track = in_track
        self.out_track = out_track
        self.smooth = int(smooth)
            
    def run(self):
        out_track_name = self.in_track._v_name + "_smooth" + str(self.smooth)
        print out_track_name
        #pdb.set_trace()
        test = tutil.checkIfNodeExists(self.out_track, out_track_name, True, False)
        if test: return
        
        for chrom in self.in_track._f_iterNodes():
            chr_name = chrom._v_name
            print chr_name
            if chr_name == "unknown": continue
            track_chr = self.in_track._f_getChild(chr_name)
            
            out_track_chr = signal_utils.smooth(track_chr[:], self.smooth, window="flat")
            
            self.out_track.createArray("/" + out_track_name, chr_name, out_track_chr)
            for name in track_chr._v_attrs._f_list():
                self.out_track.setNodeAttr("/" + "/".join([out_track_name, chr_name]), name, track_chr._v_attrs[name])
    
def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', dest='track', help='input track')
    parser.add_argument('-n', dest='subtrack', help="subtrack name. Input 'all' to process all subtracks")
    parser.add_argument('-o', dest='out_name', help='out track file')
    parser.add_argument("-s", dest="smooth", help="moving average window size. must be odd")
    args = parser.parse_args()
    
    track_file = tb.openFile(args.track)
    out_file = tb.openFile(args.out_name, 'a')
    subtrack_name = args.subtrack
  
    if subtrack_name == "all":
        for subtrack in track_file.iterNodes("/"):
            processor = trackProcess(subtrack, out_file, args.smooth)
            processor.run()
    else:
        subtrack = track_file.getNode("/" + subtrack_name)
        processor = trackProcess(subtrack, out_file, args.smooth)
        processor.run()

    out_file.flush()    
    out_file.close()
    
if __name__ == '__main__':
    main(sys.argv)
