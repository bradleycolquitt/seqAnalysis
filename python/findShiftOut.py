#!/usr/bin/env python

import sys, re
import itertools

def main(argv):
    file = open(argv[1])
    SO_ind = []
    ind = 0
    for line in file:
        if re.search("\x0e", line):
            SO_ind.append(ind)
            print line
        ind = ind + 1
    file.close()
    SO_ind_expand = []
    for i in SO_ind:
        SO_ind_expand.append([i-3, i-2, i-1, i])
    SO_ind_expand_join = itertools.chain.from_iterable(SO_ind_expand)
    ind = 0
    file = open(argv[1])
    outfile = open(argv[1] + "_trim", 'w')
    for line in file:
        if not ind in SO_ind_expand_join:
            outfile.write(line)
    file.close()
    outfile.close()
    
if __name__ == '__main__':
    main(sys.argv)