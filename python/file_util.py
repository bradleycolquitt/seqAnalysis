#!/usr/bin/env python

import sys, os, re
from subprocess import Popen

def transfer_to_dir(file_path, storage_dir):
    file_path = os.path.abspath(file_path)
    file_path_split = file_path.split("/")
    base_ind = 0
    for i in range(len(file_path_split)):
        if re.search("storage", file_path_split[i]) or re.search("user", file_path_split[i]):
            base_ind = i + 1
            
    file_base = "/".join(file_path_split[base_ind:])
    new_path = ""
    if re.search("storage", storage_dir):
        new_path = "/".join(["/media", storage_dir, file_base])
    elif re.search("home", storage_dir):
        new_path = "/".join(["/home/user", file_base])
    print "Transfer " + file_path + " to "  + new_path
    
    try:
        cmd_args = ['cp', '-r', file_path, new_path]
        p = Popen(cmd_args)
    except:
        print str(sys.exc_info()[0])
        raise
    else:
        cmd_args = ['rm', file_path]
        

    