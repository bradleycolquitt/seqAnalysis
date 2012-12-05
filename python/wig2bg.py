#!/usr/bin/env python


# Input fixedStep wig
# Output bed file with: chr start stop 0[dummy name] value[from wig] "+"[dumm strand]
# 
import re, os, shutil, time, sys
from string import *
from optparse import OptionParser
import operator

def main(argv):

	parser = OptionParser()

	parser.add_option("-i", "--input", action="store", type="string", dest="input", metavar="<str>")
	parser.add_option("-o",action="store",type="string",dest="output",metavar="<str>")
	parser.add_option("-s", action="store", type="int", dest="style")
	
	(opt, args) = parser.parse_args(argv)
	
	infile = open(opt.input, 'r');
	outfile = open(opt.output, 'w');
	
	chrom = ""
	curr_start = 0
	curr_end = 0
	step = 0
	step_pattern = re.compile("Step")
	out_line = ""
	for line in infile:
		if step_pattern.search(line):
			sline = line.split()
			chrom = sline[1].split("=")[1]
			curr_start = int(sline[2].split("=")[1])
			step = int(sline[3].split("=")[1])
			curr_end = curr_start + step - 1
			print chrom
			continue
		value = line.strip()
		if opt.style == 4 and float(value) == 0.0:
			curr_start = curr_start + step
			curr_end = curr_start + step - 1
			continue
		if opt.style == 6:
			out_line = "\t".join([chrom, str(curr_start), str(curr_end), "0", value, "+"]) + "\n"
		elif opt.style == 4:
			out_line = "\t".join([chrom, str(curr_start), str(curr_end), value]) + "\n"
		outfile.write(out_line)
		curr_start = curr_start + step
		curr_end = curr_start + step - 1
	outfile.close()
	


	


if __name__ == "__main__":
	main(sys.argv) 	
