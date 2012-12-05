#!/usr/bin/env python

import os
from subprocess import Popen

def make_tdf(wig, tdf):
    cmd_args = ['igvtools', 'tile', wig, tdf, 'mm9']
    fnull = open(os.devnull, 'w')
    Popen(cmd_args, stdout=fnull, stderr=fnull)
    fnull.close()