#!/usr/bin/env python

import sys, os
from subprocess import Popen

def main(argv):
    beds = os.listdir(argv[1])
    feature = os.path.basename(argv[2])
    out_dir = "/".join([argv[1], feature])
    if not os.path.exists(out_dir): os.mkdir(out_dir)
    
    for bed in beds:
        out = open("/".join([out_dir, bed]), 'w')
        cmd_args = ['intersectBed',
                    '-a', argv[1] + "/" + bed,
                    '-b', argv[2],
                    '-u']
        
        p = Popen(cmd_args, stdout=out)
        p.wait()
    
if __name__ == '__main__':
    main(sys.argv)