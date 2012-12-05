#! /usr/bin/env python

import sys, os
import strandBias
from multiprocessing import Process, Pool, Queue

roi_path = "/seq/lib/roi/strand/"

d3a_names = ['hmc/moe_wt_hmc.bed', 'hmc/moe_d3a_hmc.bed', 
             'mc/moe_wt_mc.bed', 'mc/moe_d3a_mc.bed',
             'in/moe_wt_in.bed', 'in/moe_d3a_in.bed']
d3a = []
for d3a_name in d3a_names:
    d3a.append("dnmt3a/" + d3a_name)

cells_names = ['hmc/omp_hmedip.bed', 'hmc/ngn_hmedip.bed', 'hmc/icam_hmedip.bed',
               'mc/omp_medip.bed', 'mc/ngn_medip.bed', 'mc/icam_medip.bed']
cells = []
for cells_name in cells_names:
    cells.append("cells/" + cells_name)

rois = os.listdir(roi_path)

def runpool(sample, roi_queue):
    pool = Pool(processes=2)
    while not roi_queue.empty():
        roi = roi_queue.get()
        print "-- " + roi
        pool.apply_async(strandBias.main, (['x', sample, roi],))
    pool.close()
    pool.join()

def main(argv):
    if argv[1] == "d3a":
        samples = d3a
    elif argv[1] == "cells":
        samples = cells
    roi_queue = Queue()
    for sample in samples:
        print sample
        for roi in rois:
            roi_queue.put(roi)
        p = Process(target=runpool, args=(sample, roi_queue))
        p.start()
        p.join()
    
if __name__ == "__main__":
    main(sys.argv)
