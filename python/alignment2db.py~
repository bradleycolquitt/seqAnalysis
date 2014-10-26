#! /usr/bin/env python 

##############################################################################################
# Parse CompletedJobInfo.xml and GenerateFASTQRunStatistics.xml within Illumina run folder
# Export
#   - library ID <SampleID>
#   - run ID <Run>
#   - instrument <Instrument>
#   - cycle number <NumCycles>
#   - type of sequencing PE/SE (Calculate from number of children of <Read>)
#   - Total number of clusters <NumberOfClustersRaw>
#   - Passed number of clusters <NumberOfClustersPF>
##############################################################################################


import sys
import xml.etree.ElementTree as ET
import MySQLdb as mdb
import pdb
import datetime

def main(argv):

    run_info = "/".join([argv[1], "CompletedJobInfo.xml"])
    fastq_info = "/".join([argv[1], "GenerateFASTQRunStatistics.xml"])
    run_file = ET.parse(run_info)
    fastq_file = ET.parse(fastq_info)
    
    ## Connect to db
    try:
        conn = mdb.connect('localhost', 'brad', 'Eu23ler1', 'sample_db')
        cur = conn.cursor()
    except mdb.Error, e:
        print "MySQLdb error %d: %s " % (e.args[0] + e.args[1])
    
    # Parse general run data
    root = run_file.getroot()
    run_node = root.find("RTARunInfo").find("Run")
    
    run_dict = dict()    
    run_dict['run_id'] = run_node.attrib['Id']
    run_dict['machine'] = run_node.find("Instrument").text
    run_dict['date'] = "20" + run_node.find("Date").text + "000000" ## Formatted for MySQL
    
    num_reads = 0
    for read in run_node.find("Reads").iter("Read"):
        num_reads += 1
        if num_reads == 1:
            run_dict['cycles'] = read.attrib['NumCycles']
    
    if num_reads == 3:
        run_dict['ended'] = "pe"
    else:
        run_dict['ended'] = "se"

    for i in run_dict.iterkeys():
        run_dict[i] = "'" + str(run_dict[i]) + "'"    
    
    # Parse FASTQ data
    root = fastq_file.getroot()
    info = dict()
    for sample in root.iter("SummarizedSampleStatistics"):
        # Populate dict
        info['idlibrary'] = sample.find("SampleID").text
        info['clusters_pf'] = sample.find("NumberOfClustersPF").text
        info['clusters_raw'] = sample.find("NumberOfClustersRaw").text
        
        # Process dict
        info['idlibrary'] = filter(lambda x: x.isdigit(), info['idlibrary'])
        for i in info.iterkeys():
            info[i] = "'" + str(info[i]) + "'"
            
        #pdb.set_trace()
        #Write to db            
        ## Check if record exists using run_id and idlibrary
        try:
            exec_string = "SELECT EXISTS(SELECT 1 FROM sequencing WHERE run_id={0} AND idlibrary={1})".format(run_dict['run_id'], info['idlibrary'])
            cur.execute(exec_string)
            res = cur.fetchall()[0]
            if res[0] == 0:
                exec_string = "INSERT INTO sequencing (" + ", ".join(run_dict.keys() + info.keys()) + ") VALUES (" + ", ".join(run_dict.values() + info.values()) + ")"
                cur.execute(exec_string)
            else:
                exec_string = "UPDATE sequencing SET " +\
                              ", ".join(["=".join(items) for items in run_dict.items()] + ["=".join(items) for items in info.items()]) +\
                              " WHERE run_id={0} AND idlibrary={1}".format(run_dict['run_id'], info['idlibrary'])
                cur.execute(exec_string)
            conn.commit()
        except mdb.OperationalError, e:
            print e
        
    #Close connection
    conn.close()

if __name__ == "__main__":
    main(sys.argv)
