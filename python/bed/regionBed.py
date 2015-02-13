#! /usr/bin/env python

import os
import shutil
import sys
import argparse
import pdb
from string import *
from math import *

"""
Takes 6 field bed
  1) chr
  2) start
  3) stop
  4) name
  5) score
  6) strand

Splits each region into windows of given size. Score is fraction of total read length
"""

def main(argv):
  parser = argparse.ArgumentParser(description="Takes five or six column file with: chr |" + 
                                   " [start] position | [end] name | score | strand \n" + 
                                   "Generates file with windows of given size around each position")
  parser.add_argument("-i", required=True, dest='input')
  parser.add_argument("-r", required=True, dest='num_regions')
  #parser.add_argument("-f", required=True, dest='flanks', help="Number of windows on each side of point")
  #parser.add_argument("-d", required=False, dest='direct', default="both", help="Which Direction from mid windows should be made")
  args = parser.parse_args()
  infile = open(args.input, 'r')
  #outfile = open(args.input + "_W" + args.window + "F" + args.flanks + "_" + args.direct, 'w')
  outfile = open(args.input + "_R" + args.num_regions, 'w')
  regions = atoi(args.num_regions)
  #flanks = atoi(args.flanks)
  #direct = args.direct
  genome_sizes = open('/seq/lib/mouse.mm9.genome', 'r')

  d = {}
  for line in genome_sizes:
      line = line.strip()
      sline = line.split()
      if len(sline) > 0:
          d[sline[0]]=sline[1]
  strand = "+"
  name = "None"
  
  for line in infile:
    #pdb.set_trace()
    sline = line.strip().split();
    if len(sline) == 4:
      name = sline[3]
      strand = "+"
    if len(sline) == 5:
      name = sline[2]
      strand = sline[4]
    elif len(sline) == 6:
      name = sline[3]
      strand = sline[5]
    if strand == "+":
      start = int(sline[1])
      end = int(sline[2])
    else:
      start = int(sline[2])
      end = int(sline[1])
    bed_length = int(sline[2]) - int(sline[1])
    window = bed_length / regions
    if window == 1: continue
    r = range(regions)
  #    if strand == "+":
  ##      if direct == "up" or direct=="both":
  #      if len(sline) == 5:
  #       start = atoi(sline[1]) - flanks * window + 1
  #      elif len(sline) == 4 or len(sline) == 6:
  #        start = ((int(sline[1]) + int(sline[2])) / 2) - flanks * window + 1
  #      #elif direct=="down":
  #      if len(sline) == 5:
  #       start = atoi(sline[1])
  #      elif len(sline) == 4 or len(sline) == 6:
  #        start = ((int(sline[1]) + int(sline[2])) / 2)
  #    elif strand == "-":
  #      if direct == "up" or direct=="both":
  #        if len(sline) == 5:
  #          end = atoi(sline[1]) + flanks * window 
  #        elif len(sline) == 6:
  #          end = ((int(sline[1]) + int(sline[2])) / 2) + flanks * window
  #      elif direct=="down":
  #        if len(sline) == 5:
  #         end = atoi(sline[1])
  #        elif len(sline) == 6:
  #          end = int(sline[2]) - ((int(sline[1]) + int(sline[2])) / 2)
      
  #    if direct == "both":
  #      r = xrange(2 * flanks)
  #    else:
  #      r = xrange(flanks)
    out = ""
    for index in r:
      out = (sline[0] + "\t" + str(start) + "\t" + str(end) + "\t" + \
                name + "\t" + str(index) + "\t" + strand + "\n")
      if (start and end >= 0) and (start and end <= atoi(d[sline[0]])):
          outfile.write(out)
      if strand == "+":
          start = start + window
          end = end + window
      else:
          start = start - window        
          end = end - window
  outfile.close()

if __name__ == "__main__":
    main(sys.argv)
