#! /bin/awk

#@1 input fastq
# @2 length of output

awk {
    if NR % 2 == 0 {
	    print substr($1, 1, )
}
}
