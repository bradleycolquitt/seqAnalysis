#!/usr/bin/env python

import re, os, shutil, time, sys
from math import *
from string import *
from optparse import OptionParser
import operator

#import BED_mod
#import GenomeData

def main(argv):

	parser = OptionParser()
	parser.add_option("-i", "--infile", action="store", type="string", dest="infile", metavar="<str>")
	#parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", metavar="<file>")
	parser.add_option("-t", action="store_true", dest="title", metavar="<boolean>",default=False)
	#parser.add_option("-g", action="store",type="string",dest="genome",metavar="<str>")
	
	(opt, args) = parser.parse_args(argv)
	infile = open(opt.infile, 'r');
      	outfile_name = opt.infile + "_groom" 
	outputfile = open(outfile_name, 'w')	
	d={}
	genome = open("/seq/lib/mouse.mm9.genome", 'r')
	#genome = open(opt.genome,'r')
	for line in genome:
		line = line.strip()
		sline = line.split()
		if len(sline)>0:
			d[sline[0]] = sline[1]

	if opt.title:
		out = "track name=\"" + str(opt.infile) + "\"\t" + "description=\"" + str(opt.outfile) + "\"\n"	
		outputfile.write(out)

	for line in infile:
		line = line.strip()
		sline = line.split()
		if atoi(sline[2]) <= atoi(d[sline[0]]):		
			out = "\t".join(sline) + "\n"
			outputfile.write(out)
	outputfile.close()

if __name__ == "__main__":
	main(sys.argv)
