	#!/usr/bin/env python

import re, os, shutil, time, sys
from string import *
from optparse import OptionParser
import operator

 

def main(argv):

	parser = OptionParser()

	parser.add_option("-i", "--input", action="store", type="string", dest="input", metavar="<str>")
	parser.add_option("-l",action="store",type="int",dest="select",metavar="<int>")
	
	(opt, args) = parser.parse_args(argv)
	
	infile = open(opt.input, 'r');
	select = opt.select
	sum = 0;

	for line in infile:
		sum += 1;
		
	
	print sum;

	


if __name__ == "__main__":
	main(sys.argv) 	
