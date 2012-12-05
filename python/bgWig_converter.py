#!/usr/bin/env python

import re, os, sys
import argparse
from string import *

def bg2wig(infile, outfile):
	chrom_store = ""
	for line in infile:	
		line = line.strip()
		sline = line.split()
		chrom = sline[0]
		start = atoi(sline[1])
		end = atoi(sline[2])
		val = sline[3]
		step = str(end - start + 1)
		
		if chrom == chrom_store:
			out = val + "\n"
			outfile.write(out)
		else:	
			out = " ".join(["fixedStep", "chrom=" + chrom, "start=" + str(start), "step=" + str(step)]) + "\n"
			outfile.write(out)
			chrom_store = chrom

def wig2bg(infile, outfile):
	header_pattern = re.compile("Step")
	step_type = ""
	chrom = ""
	curr_pos = 0
	next_pos = 0
	step = 0
	for line in infile:
		line = line.split()
		if re.search(header_pattern, line[0]):
			step_type = line[0]
			fields = []
			for field in line[1:]:
				fields.append(field.split("=")[1])
			chrom = fields[0]
			if step_type == "fixedStep":
				curr_pos = atoi(fields[1])
				step = atoi(fields[2])
		else:
			if step_type == "fixedStep":
				val = line[0]
				next_pos = curr_pos + step - 1
				out = "\t".join([chrom, str(curr_pos), str(next_pos), val]) + "\n"
				outfile.write(out)
				curr_pos = next_pos + 1


def main(argv):
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", dest="input")
	args = parser.parse_args()
	
	input_prefix = args.input.split(".")[0]
	infile = open(args.input, 'r')

	if re.search("wig", args.input):
		outfile = open(input_prefix + ".bed", 'w')
		wig2bg(infile, outfile)
	elif re.search("bed", args.input):
		outfile = open(input_prefix + ".wig", 'w')
		bg2wig(infile, outfile)
	infile.close()
	outfile.close()
if __name__ == "__main__":
	main(sys.argv)
