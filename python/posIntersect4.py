#!/usr/bin/env python

#Import two bed files
#	1.summary BED(chr,start,stop,fold_change)
#	2.gene BED(chr,start,stop,id)
#Set up export for each chromosome
#For each line in each chromosome block of summary BED,
#	 find interval in gene bed such that:
#		summary.start>gene.start-1000
#		AND
#		summary.stop<gene.start+1000
#	 write to export (id,value)

import re, os, shutil, time, sys
from math import *
from string import *
from optparse import OptionParser
import operator

import BED_mod
import GenomeData
import lineCount

"""
Count number of lines in gene file
For each test
        bisect gene file
        if intersect: write
        elif greater: bisect above
        elif lesser: bisect below
"""


def main(argv):

	parser = OptionParser()

	#parser.add_option("-s", "--species", action="store", type="string", dest="species", help="mm8, hg18, background, etc", metavar="<str>")
	parser.add_option("-t", "--test", action="store", type="string", dest="test", metavar="<str>")
	parser.add_option("-g", "--gene", action="store", type="string", dest="gene", metavar="<str>")
	parser.add_option("-i", "--overlap", action="store", type="int", dest="overlap", metavar="<int>",default=0)
	parser.add_option("-o", "--outFile", action="store", type="string", dest="outFile", metavar="<file>")
	parser.add_option("-k", "--kind", action="store", type="string", dest="kind",metavar="<str>")
	parser.add_option("-r", "--reads", action="store_true", dest="reads",metavar="<boolean>",default=False)
	#parser.add_option("-f", "--fraction", action="store", type="float", dest="fraction", metavar="<float>",default=0.0)

	
	(opt, args) = parser.parse_args(argv)
	
	overlap=opt.overlap;
	#fraction=opt.fraction;
	#if opt.species in GenomeData.species_chroms.keys():
	#	print "Species: ", opt.species;
	gene = open(opt.gene,'r')
	test = open(opt.test,'r')
        gn = ""
	tn = ""
	gline,tline = gene.readline(),test.readline();
	sgline,stline = gline.split(),tline.split();
        bed_length = [len(sgline),len(stline)]
        
        types = [1,2];
        
        for index in xrange(len(bed_length)):
                if bed_length[index] == 4: types[index] = "BED_GRAPH"    
                elif bed_length[index] == 6: types[index] ="BED6"
        if opt.reads: types[0]="BED7"
        
	if type(sgline[3]) == type("x"):
		gn = "name"
	if type(stline[3]) == type("x"):
		tn = "name"
		

        
	print "Accessing test BED...";
	testBED = BED_mod.BED(species="mm9",file=opt.test,bed_type=types[1],n=tn);
	
	print "Accessing gene BED...";
	geneBED = BED_mod.BED(species="mm9",file=opt.gene,bed_type=types[0]);
	
	outputfile = open(opt.outFile, 'w');
	
	sums = 0;
	
	for chrom in testBED.keys():

                count = len(geneBED.bed_vals[chrom]);
                if count == 0: continue;

		print "Test chromosome: " + chrom;
		gstarts = geneBED.getStarts(chrom)
		print gstarts
		sys.exit()
		if opt.kind == "in":
			for t in testBED.bed_vals[chrom]:

                                
                                bottom = 0;
                                top = count-1;
                                switch = 1;
                                bisect = int(round(count/2)-1)
                                
                                while switch:
                                        
                                        g = geneBED.bed_vals[chrom][bisect];
                                       
                                        if conditions(t.start,t.end,g.start,g.end,overlap):
                                                
                                                sums += 1;
						g.reads += 1;
                                                switch = 0;
                                                
                                        elif t.start > g.start:
                        
                                                bottom = bisect
                                                bisect_t = int(round((top-bottom+1)/2) + bottom - 1)
						
                                                if bisect_t == bisect: switch = 0;
                                                else: bisect = bisect_t
                                        elif t.end < g.end:
                                                
                                                top = bisect
                                                bisect_t = int(round((top-bottom+1)/2) + bottom - 1)
                                                
						if bisect_t == 0: switch = 0;
						else: bisect = bisect_t
			for chrom in geneBED.keys():
                                for val in geneBED.bed_vals[chrom]:
                                        pos = val.getCoord() + "\n"
                                        outputfile.write(pos);		
			
		elif opt.kind == "out":
			for t in testBED.bed_vals[chrom]:
                                
                                bottom = 0;
                                top = count-1;
                                switch = 1;
                                bisect = int(round(count/2)-1)
                                out = ""
                                while switch:
                                        
                                        g = geneBED.bed_vals[chrom][bisect];
                                       
                                        if conditions(t.start,t.end,g.start,g.end,overlap):
                                                
                                                #sums += 1;
						#g.reads += 1;
                                                switch = 0;
                                                
                                        elif t.start > g.start:
                        
                                                bottom = bisect
                                                bisect_t = int(round((top-bottom+1)/2) + bottom - 1)
						
                                                if bisect_t == bisect:
                                                        switch = 0;
                                                        sums += 1;
                                                        out = t.getCoord() + "\n"
                                                else: bisect = bisect_t
                                        elif t.end < g.end:
                                                
                                                top = bisect
                                                bisect_t = int(round((top-bottom+1)/2) + bottom - 1)
                                                
						if bisect_t == 0:
                                                        switch = 0;
                                                        sums += 1;
                                                        out = t.getCoord() + "\n"
						else: bisect = bisect_t
				outputfile.write(out);
						
				
	print "Number of hits: " + str(sums)
	outputfile.close();

def index(a, x):
    'Locate the leftmost value exactly equal to x'
    i = bisect_left(a, x)
    if i != len(a) and a[i] == x:
        return i
    raise ValueError


def conditions(ts,te,gs,ge,o):
	a = 0;
	if ts<=gs and te>=ge:
		a = 1;	
	elif ts>=gs and te<=ge: 
		a = 1;
		
	elif ts<=gs and te>=gs+o:
		a = 1;

	elif ts<=ge-o and te>=ge: 
		a = 1;
	
	return a; 	
	
						

        
if __name__ == "__main__":
	main(sys.argv)
