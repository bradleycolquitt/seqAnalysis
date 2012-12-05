#!/usr/bin/env python

import sys, os, re
import itertools

def main(argv):
    files = os.listdir(argv[1])
    term = re.compile("fastq$")
    files = filter(term.search, files)
    print files
    for file_name in files:
        file = open(file_name)
        SO_ind = []
        ind = 0
        terms = []
        for i in [14, 17, 18]:
            terms.append(chr(i))
        print terms
        #sys.exit()
        for line in file:
            for term in terms:
                if re.search(term, line):
                    SO_ind.append(ind)
                    print line
            ind = ind + 1
        file.close()
        
        SO_ind_expand = []
        for i in SO_ind:
            SO_ind_expand.extend([i-3, i-2, i-1, i])
    #    SO_ind_expand_join = itertools.chain.from_iterable(SO_ind_expand)
        print SO_ind_expand
        #sys.exit()
        ind = 0
        file = open(file_name)
        outfile = open(file_name + "_trim", 'w')
        for line in file:
            if not ind in SO_ind_expand:
                outfile.write(line)
            ind = ind + 1    
        file.close()
        outfile.close()
    
if __name__ == '__main__':
    main(sys.argv)