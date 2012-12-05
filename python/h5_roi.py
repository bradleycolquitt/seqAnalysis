#!/usr/bin/env python
from string import *
from tables import *
import numpy
import re

class ROI(IsDescription):
    name = StringCol(20)
    window = UInt8Col()
    number = UInt8Col()
    flank = UInt8Col()
    class ROI_chrs(IsDescription):
            _v_pos = 1
            chr = StringCol(2)
            start = UInt32Col()
            end = UInt32Col()
            name = StringCol(16)
            pos = Int16Col()
            strand = StringCol(1)
    
def loadBed(bed_name, leaf):
    bed = open(bed_name)
    row = leaf.row
    name =  bed_name.split("_")
    w = re.compile("W")
    info = filter(w.search, name)[0]
    info = info.split("[A-Z]")
    row = leaf.ROI_chrs.row
    for line in bed:
        line = line.strip().split()
        row['ROI_chrs/chr'] = line[0].split('chr')[1]
        row['ROI_chrs/start'] = line[1]
        row['ROI_chrs/end'] = line[2]
        row['ROI_chrs/name'] = line[3]
        row['ROI_chrs/pos'] = line[4]
        row['ROI_chrs/strand'] = line[5]
        row.append()
    leaf.flush()