#!/usr/bin/env python

def checkIfNodeExists(out, name):
    nodes = out.listNodes("/")
    try:
        out.getNode("/" + name)
    except tb.NoSuchNodeError:
        pass
    else:
        dec = raw_input("Node exists. Overwrite [y/n]? ")
        if dec == "n": return True
        elif dec == "y":
            out.removeNode("/" + name, recursive=True)
    out.createGroup("/", name)
    return False