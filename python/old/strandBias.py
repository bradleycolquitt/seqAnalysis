#!/usr/bin/env python

import re, os, shutil, time, sys, getpass
from string import *
from subprocess import Popen
from numpy import array, zeros

DATA_PATH = "/media/storage2/data/working"
ROI_PATH = "/seq/lib/roi/strand"
OUT_PATH = "/media/storage2/analysis/strand"

class bed_sample:
  def __init__(self, sample):
    self.name = sample
    self.root = split(os.path.basename(sample), ".bed")[0]
    self.path = "/".join(split(sample, "/")[0:-1])
    self.plus_name = DATA_PATH + "/" + self.path + "/" + self.root + "_plus"
    self.minus_name = DATA_PATH + "/" + self.path + "/" + self.root + "_minus"
    
class roi:
  def countRoiLines(self):
    nlines = 0
    nfields = 0
    roi_file = open(self.path, 'r')
    for line in roi_file:
      nfields = len(line.strip().split())
      nlines = nlines + 1
    return (nlines, nfields)
    
  def __init__(self, roi):
    self.name=roi
    self.root=split(roi, ".bed")[0]
    self.path = ROI_PATH + "/" + roi
    self.nlines=self.countRoiLines()[0]
    self.nfields=self.countRoiLines()[1]
    
class strand:
  def __init__(self, sample_name, roi_name):
    self.sample = bed_sample(sample_name)
    self.roi = roi(roi_name)
    self.out_path = ""
    self.plus_out = ""
    self.minus_out = ""
    
  def split(self):
    bed = open(DATA_PATH + "/" + self.sample.name, 'r')
    if not os.path.exists(self.sample.plus_name):
      bed_plus = open(self.sample.plus_name, 'w')
      bed_minus = open(self.sample.minus_name, 'w')
      print "Splitting plus and minus reads"
      for line in bed:
        if re.search("\+", line):
          bed_plus.write(line)
        else:
          bed_minus.write(line)
      bed_plus.close()
      bed_minus.close()

  def intersectRoi(self):
    self.out_path = OUT_PATH + "/" + self.sample.root
    if not os.path.exists(self.out_path):
      os.mkdir(self.out_path)
    self.plus_out = self.out_path + "/" + self.roi.root + "_plus"
    self.minus_out = self.out_path + "/" + self.roi.root + "_minus"
    if os.path.exists(self.plus_out):
      return
    print "Intersecting reads with ROI"
    plus_out = open(self.plus_out, 'w')
    minus_out = open(self.minus_out, 'w')
    cmd_args_plus = ['intersectBed', '-a', ROI_PATH + "/" + self.roi.name, 
                     '-b', self.sample.plus_name, '-wa', '-c']
    cmd_args_minus = ['intersectBed', '-a', ROI_PATH + "/" + self.roi.name, 
                      '-b', self.sample.minus_name, '-wa', '-c']
  
    cmds = [cmd_args_plus, cmd_args_minus]
    outs = [plus_out, minus_out]
    proc = []
    for i in xrange(len(cmds)):
      proc.append(Popen(cmds[i], stdout=outs[i]))
    for p in proc:
      p.wait()
    plus_out.close()
    minus_out.close()
  
  def computeScore(self):
    print "Computing strand scores"
    score_position = 0
    strand_position = 0
    if self.roi.nfields == 5: 
      score_position = 5
      strand_position = 3
    elif self.roi.nfields == 6: 
      score_position = 6
      strand_position = 5 
    vals = zeros((self.roi.nlines,3))
    index = 0
    for line in open(self.roi.path, 'r'):
      line = line.split()
      vals[index, 0] = strandToNum(line[strand_position])
      index = index + 1
    index = 0
    for line in open(self.plus_out, 'r'):
      line = line.split()
      vals[index,1] = atoi(line[score_position])
      index = index + 1
    index = 0
    for line in open(self.minus_out, 'r'):
      line = line.split()
      vals[index,2] = atoi(line[score_position])
      index = index + 1
    #return vals
    scores = list()
    index = 0
    #return vals
    for row in vals:
      try:
        scores.append(compareVals(row))
        index = index + 1
      except:
        print "compareVals failed:", sys.exc_info()[0]
        print "row value:", row
        print "return value:", scores[index]
        raise
    return scores

def strandToNum(val):
  out = 2
  if val == "+":
    out = 1
  elif val == "-":
    out = 0
  return out
    
def compareVals(row):
  score = 0
  denom = float(row[1] + row[2])
  if denom > 0:
    if row[0] == 1:
      score = float(row[1] - row[2]) / denom
    elif row[0] == 0:
      score = float(row[2] - row[1]) / denom
  return [str(denom), str(score)]

def main(argv):
  s = strand(argv[1], argv[2])
  #print "strand object created"
  s.split()
  #print "split"
  s.intersectRoi()
  #print "intersected"
  score_path = s.out_path + "/" + s.roi.root + "_scores.txt"
  if os.path.exists(score_path):
    #print "returning" + score_path
    return
  score_file = open(score_path, 'w')  
  try:
    scores = s.computeScore()
  except:
    print "Scores not computed: ", sys.exc_info()[0]
    print "-- file =", s.roi.root
    shutil.rmtree(score_path)
    raise
  #print "scores computed"
  #print "before open"
  roi = open(s.roi.path, 'r')
  #print "after open"
  write_position = 0
  if s.roi.nfields == 5: write_position = 5
  elif s.roi.nfields == 6: write_position = 6
  #print write_position
  index = 0
  for line in roi:
    line = line.strip().split()
    out = "\t".join(line[:write_position]) + "\t" + "\t".join(scores[index]) + "\n"
    score_file.write(out)  
    index = index + 1
  #print "scores written"
  score_file.close()
  
if __name__ == "__main__":
  main(sys.argv)
