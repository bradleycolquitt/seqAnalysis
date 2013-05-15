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
    def __init__(self, in_track, out_name, dampen):
        self.in_track = in_track
        self.in_track_name = in_track._v_name
        self.out_name = self.in_track_name + "_dampen" + dampen + ".bed"
        self.out_file = open(self.out_name, 'w')
        self.dampen = int(dampen)
            
    def run(self):
        # various machinations to keep last name component last
        #track_name_split = self.in_track._v_name.split("_")
        #out_track_name = "_".join(["_".join(track_name_split[0:-1]), "dampen" + str(self.dampen), track_name_split[-1]])
        #print out_track_name
        #pdb.set_trace()
        #test = tutil.checkIfNodeExists(self.out_track, out_track_name, True, False)
        #if test: return
        L = list(self.in_track._f_iterNodes())
        cutoff = [0]
        #combined = np.zeros(1)
        #pdb.set_trace()
        for chrom in L: 
            chr_name = chrom._v_name
            print chr_name
            track_chr = self.in_track._f_getChild(chr_name)
#            combined = np.concatenate(combined, track_chr)
        
            cutoff.append((self.dampen * np.std(track_chr[:])))
        cutoff = max(cutoff)    
        #pdb.set_trace()

        for chrom in iter(L):
            chr_name = chrom._v_name
            track_chr = self.in_track._f_getChild(chr_name)
            mask = ~(np.abs(track_chr[:] - np.mean(track_chr[:])) < cutoff)
            self.coord_to_bed(mask, chr_name, track_chr.attrs.Stop, track_chr.attrs.Resolution, self.out_file)
            #out_track_chr = track_chr[:]
            #out_track_chr[mask] = cutoff        
            #self.out_track.createArray("/" + out_track_name, chr_name, out_track_chr)
            #for name in track_chr._v_attrs._f_list():
            #    self.out_track.setNodeAttr("/" + "/".join([out_track_name, chr_name]), name, track_chr._v_attrs[name])
    
    def coord_to_bed(self, mask, chrom, stop, res, out):
        start = 0
        end = 0
        #pdb.set_trace()
        coord = np.r_[1:(stop + 1)]
        coord = coord[mask]
        for pos in coord:
            start = int(((pos - 1) * res) + 1)
            end = int(start + (res - 1))
            out.write("\t".join([chrom, str(start), str(end), str(pos), "0", "+"]) + "\n")

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', dest='track', help='input track')
    parser.add_argument('-n', dest='subtrack', help="subtrack name. Input 'all' to process all subtracks")
    parser.add_argument('-o', dest='out_name', help='out bed file')
    parser.add_argument("-d", dest="dampen", help="number of standard deviations away from mean before dampening")
    args = parser.parse_args()
    
    track_file = tb.openFile(args.track)
    #out_file = open(args.out_name, 'w')
    subtrack_name = args.subtrack
  
    if subtrack_name == "all":
        for subtrack in track_file.iterNodes("/"):
#            processor = trackProcess(subtrack, out_file, args.dampen)
            processor = trackProcess(subtrack, args.out_name, args.dampen)
            processor.run()
    else:
        subtrack = track_file.getNode("/" + subtrack_name)
#        processor = trackProcess(subtrack, out_file, args.dampen)
        processor = trackProcess(subtrack, args.out_name, args.dampen)
        processor.run()

    out_file.flush()    
    out_file.close()
    
if __name__ == '__main__':
    main(sys.argv)
