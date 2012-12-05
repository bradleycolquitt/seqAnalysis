#!/usr/bin/env python
import os
from Bio import SeqIO
from multiprocessing import Pool

count1 = 0
count2 = 0

def worker(args):
    input = args[0]
    read1 = args[1]
    read2 = args[2]
    file = args[3]
    count = 0
    print file
    indata = open("/".join([input,file]))
    file_split = file.split("_")
    if int(file_split[2]) == 1:
        out = read1 
    elif int(file_split[2]) == 3:
        out = read2
    elif int(file_split[2]) == 2:
        return
    for line in indata:
        count = count + 1
        line = line.strip().split()
        head = "%s:%s:%s:%s:%s:%s#%s/%s PF=%s" % tuple(line[0:8] + [line[10]])
        read = line[8].replace(".", "N")
        qual = line[9]
        ok = True
        for c in qual:
            if ord(c) < 64:
                ok = False
        if ok:
            out.write("\n".join(["@" + head, read, "+", qual]) + "\n")
    indata.close()
    return count

def convertToFastq(dirs):
    input = dirs[0]
    output = dirs[1]
    if not os.path.exists(input):
        return
    if not os.path.exists(output): os.mkdir(output)
    else: return
                     
    files = os.listdir(input)
    read1 = open("/".join([output, "1.fastq"]), 'w')
    read2 = open("/".join([output, "2.fastq"]), 'w')
    count1 = 0
    count2 = 0
    args = [(input, read1, read2, file) for file in files]
    #print args
    #pool = Pool(10)
    for arg in args:    
        #pool.apply_async(worker, (arg,))
        count = worker(arg)
        if arg[3].split("_")[2] == "1":
            count1 = count1 + count
        elif arg[3].split("_")[2] == "3":
            count2 = count2 + count
    #pool.close()
    #pool.join()
    
    #for file in files:
    #    count = 0
    #    print file
    #    indata = open("/".join([input,file]))
    #    file_split = file.split("_")
    #    if int(file_split[2]) == 1:
    #        out = read1 
    #    elif int(file_split[2]) == 3:
    #        out = read2
    #    elif int(file_split[2]) == 2:
    #        continue
    #    for line in indata:
    #        count = count + 1
    #        line = line.strip().split()
    #        head = "%s:%s:%s:%s:%s:%s#%s/%s PF=%s" % tuple(line[0:8] + [line[10]])
    #        read = line[8].replace(".", "N")
    #        qual = line[9]
    #        out.write("\n".join(["@" + head, read, "+" + head, qual]) + "\n")
    #    if int(file_split[2]) == 1:
    #        count1 = count1 + count
    #    elif int(file_split[2]) == 3:
    #        count2 = count2 + count
    #    indata.close()
    count1_file = open("/".join([output, "1.fastq_count"]), 'w')
    count2_file = open("/".join([output, "2.fastq_count"]), 'w')
    count1_file.write(str(count1))
    count2_file.write(str(count2))
    read1.close()
    read2.close()
    count1_file.close()
    count2_file.close()
    
    
#def illumina2sanger(input, output):

    