#!/usr/bin/env python

# Convert back to BAM
    # Sort by query name
    #sort_bam = bam_prefix + "_qsort"
    #if not os.path.exists(sort_bam):
    #    cmd_args = ['java', '-Xmx6g', '-jar', '/seq/picard/SortSam.jar',
    #                "=".join(["INPUT", bam]), "=".join(["OUTPUT", sort_bam]),
    #                "=".join(["SORT_ORDER", "queryname"])]
    #    #cmd_args = ['samtools', 'sort', '-n', bam, sort_bam]
    #    p = Popen(cmd_args)
    #    p.wait()
    #
    ## bedtools bamToBed to convert to Bed
    #
    ##bed_name = bam_prefix + ".bed"
    ##bed_file = open(bed_name, 'w')
    #bam_frag_name = bam_prefix + "_frag.bam"
    #bam_frag_file = open(bam_frag_name, 'w')
    ##cmd_args1 = ['bamToBed', '-i', sort_bam + ".bam", '-bedpe']
    #cmd_args1 = ['bamToBed', '-i', sort_bam + ".bam"]
    ##cmd_args2 = ['cut', '-f', '1-2,6-9']
    #cmd_args3 = ['bedToBam', '-i', 'stdin', '-g', '/seq/lib/mouse.mm9.genome_norand']
    #tmp_name = sort_bam + "tmp"
    #try:
    #    print "Extending bed"
    #    #p1 = Popen(cmd_args1, stdout=PIPE)
    #    tmp2 = open(sort_bam + "_tmp2", 'w')
    #    p1 = Popen(cmd_args1, stdout=tmp2)
    #    p1.wait()
    #    #p2 = Popen(cmd_args2, stdin=p1.stdout, stdout=PIPE)
    #    tmp = open(tmp_name, 'w')
    #    stitchBed(p1.stdout, tmp)
    #    tmp.close()
    #    tmp = open(tmp_name, 'r')
    #    p3 = Popen(cmd_args3, stdin=tmp, stdout=bam_frag_file)
    #    p3.wait()
    #except AssertionError as detail:
    #    print "Assertion error ", detail
    #    return
    #except:
    #    print "General fail"
    #    raise
    #else:
    #    print "Extension completed successfully"
    #    os.remove(tmp_name)
    #finally:    
    #    bam_frag_file.close()
    #
    ## Sort by coordinate
    #try:
    #    print "Starting coordinate sort"
    #    pysam.sort(bam_frag_name, bam_prefix + "_frag_sort")
    #except SamtoolsError:
    #    print "Failed coordinate sort"
    #    return
    #else:
    #    print "Completed coordinate sort"
    #    os.remove(bam_frag_name)
    #    
    ## Replace header
    #if replace_header:
    #    cmd_args1 = ['samtools', 'view', '-h', bam]
    #    cmd_args2 = ['samtools', 'reheader', '-', bam_prefix + "_frag_sort.bam"]
    #    tmp = open(bam_prefix + "_tmp", 'w')
    #    try:
    #        p1 = Popen(cmd_args1, stdout=PIPE)
    #        p2 = Popen(cmd_args2, stdin=p1.stdout, stdout=tmp)
    #        p2.wait()
    #        tmp.close()
    #    except:
    #        print "Failed reheader"
    #        os.remove(tmp)
    #        return
    #    else:
    #        #os.remove(bam)
    #        os.rename(bam_prefix + "_tmp", bam_prefix + "_frag_sort.bam")
    #pysam.index(bam_prefix + "_frag_sort.bam")