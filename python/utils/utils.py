#!/usr/bin/env python

from subprocess import Popen

def make_tdf(wig, tdf):
    cmd_args = ['igvtools', 'tile', wig, tdf, 'mm9']
    Popen(cmd_args)