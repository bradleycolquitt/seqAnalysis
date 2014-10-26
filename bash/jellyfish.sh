#! /bin/bash

jellyfish count -m 50 -s 10M -t 6 -o jf/"$1"_counts.jf -C $1
