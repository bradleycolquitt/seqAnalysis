#! /bin/bash

# Simple wrapper for picard EstimateLibraryComplexity
java -jar /seq/picard_current/EstimateLibraryComplexity.jar INPUT=$1 OUTPUT=$1_complexity.txt
