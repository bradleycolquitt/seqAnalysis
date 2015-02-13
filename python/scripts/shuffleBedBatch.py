#!/usr/bin/env python

import sys, os
from subprocess import Popen

def main(argv):
    
    out_dir = argv[1] + "rand"
    if not os.path.exists(out_dir): os.mkdir(out_dir)
    for i in range(int(argv[2])):
        out = open(out_dir + "/rand" + str(i), 'w')
        cmd_args = ['shuffleBed',
                    '-i', argv[1],
                    '-g', '/seq/lib/mouse.mm9.genome']
        
        p = Popen(cmd_args, stdout=out)
        p.wait()
    
if __name__ == '__main__':
    main(sys.argv)