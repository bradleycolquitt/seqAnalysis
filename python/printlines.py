	#!/usr/bin/env python

import re, os, shutil, time, sys
from string import *
from optparse import OptionParser
import operator

 

def main(argv):

	parser = OptionParser()

	parser.add_option("-i", "--input", action="store", type="string", dest="input", metavar="<str>")
	parser.add_option("-l",action="store",type="str",dest="select",metavar="<str>")
	
	(opt, args) = parser.parse_args(argv)
	
	infile = open(opt.input, 'r');
	select = opt.select
	ss = select.split(',');
	sum = 1;
	int_sum = 0
	
	for line in infile:
		if sum < atoi(ss[0]):
			sum += 1;
		else: 
			
			if int_sum <= atoi(ss[1]):
				print line.strip();
				int_sum+=1;
			else: break; 

	


	


if __name__ == "__main__":
	main(sys.argv) 	
