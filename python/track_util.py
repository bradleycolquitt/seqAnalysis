#!/usr/bin/env python

import pdb
import tables as tb
import numpy as np

def checkIfNodeExists(out, name, create, accept_existence):
    #pdb.set_trace()
    nodes = out.listNodes("/")
    try:
        out.getNode("/" + name)
    except tb.NoSuchNodeError:
        if create:
            out.createGroup("/", name)
            return True
        else: pass
    else:
        if not accept_existence:
            dec = raw_input("Node exists. Overwrite [y/n]? ")
            if dec == "n": return True
            elif dec == "y":
                out.removeNode("/" + name, recursive=True)
        else:
            return True
    #out.createGroup("/", name)
    return False

def nodeExists(out, name, remove):
    nodes = out.listNodes("/")
    try:
        out.getNode("/" + name)
    except tb.NoSuchNodeError:
        return False
    else:
        return True

    
def setTrackAttributes(file, node, start, stop, name, resolution):
    file.setNodeAttr("/" + node, "AStart", np.array([round(start * resolution),], dtype="int32"))
    file.setNodeAttr("/" + node, "AStop", np.array([round(stop * resolution)], dtype="int32"))
    file.setNodeAttr("/" + node, "Name", name)
    file.setNodeAttr("/" + node, "Resolution", np.array(resolution, dtype="int32"))
    file.setNodeAttr("/" + node, "Start", np.array(start, dtype="int32"))
    file.setNodeAttr("/" + node, "Stop", np.array(stop, dtype="int32"))
    
def writeTrack(file, node, name, values, resolution):
    ## Check this
    template_h5 = tb.openFile("/media/storage2/data/h5/nuc_01234_25.trk")
    template_chr = template_h5.getNode("/omp_nuc_0123", name)
    if checkIfNodeExists(file, node, True, True):
        if nodeExists(file, "/".join([node, name]), True):
            file.removeNode("/" + "/".join([node, name]))
        file.createArray("/" + node, name, values)
        for attr_name in template_chr._v_attrs._f_list():
            file.setNodeAttr("/" + "/".join([node, name]), attr_name, template_chr._v_attrs[attr_name])
#        setTrackAttributes(file, node + "/" + name, 0, len(values), name, resolution)
        
        