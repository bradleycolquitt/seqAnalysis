#!/usr/bin/env python

import sys
import os
import re
import argparse
import pysam
import pdb
import corrGenome
import tables as tb
from subprocess import Popen

def wig2Track(wig, track, window):
    cmd_args = ['LoadData', '-i', wig,
                '-o', track,
                '-t', os.path.basename(wig).split(".wig")[0],
                '-n', window,
                '-g', '/media/storage2/genomedata/chromosomes.trk']
    p = Popen(cmd_args)
    p.wait()
      

        
def main(argv):
    
    parser = argparse.ArgumentParser()
    parser.add_argument(dest="trackfile")
    parser.add_argument('-o', dest='outfile')
    parser.add_argument('-f', dest='function')
    parser.add_argument('-w', dest='window')
    args = parser.parse_args()
    
    
    tfile = tb.openFile(args.trackfile)
    samples = [s._v_name for s in tfile.iterNodes("/")]
    samples_prefix = [s.split("_plus")[0] for s in samples if re.search("plus", s)]
    
    # Get track resolution
    res = str(tfile.getNode("/" + samples[0]).chr1.getAttr("Resolution")[0])
    #pdb.set_trace()
    
    
    for sample in samples_prefix:
        print sample
        wig_name = "/".join(["/media/storage2/analysis/genomecor/strand",
                             args.function,
                             "_".join([sample, "plus", sample, "minus", "W" + args.window, "S" + res]) + ".wig"])
        
        if not os.path.exists(wig_name):
             
            # Compute strand differences
            cg = corrGenome.corr_genome(h5=args.trackfile,
                        track1=sample + "_plus", track2=sample + "_minus",
                        data_type="strand", method=args.function,
                        window=args.window, step=res)
            corrGenome.compute_by_chrom(cg)
        
        #Load WIG to H5
        wig_name = "/".join(["/media/storage2/analysis/genomecor/strand", args.function, "_".join([sample, "plus", sample, "minus", "W" + args.window, "S" + res]) + ".wig"])
        wig2Track(wig_name, args.outfile, res)
    
    # LoadTrack
    #Wig to Track
    ##for bam in split_bam:
    ##    print "WIG to Track..."
     #   print bam
    #    wig2Track(bam.split(".bam")[0] + ".wig", args.outfile, args.window)
        
if __name__ == '__main__':
    main(sys.argv)
