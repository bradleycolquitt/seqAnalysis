#!/usr/bin/env python

import sys
from math import *
from string import *

def f(x,y): return pow(atoi(x),atoi(y))

def main(argv):
    
    #input list of numbers and return their squares
    #using list comprehensions

    arguments = argv[1:]
    print  reduce(f, arguments)


if __name__ == "__main__":
    main(sys.argv)
