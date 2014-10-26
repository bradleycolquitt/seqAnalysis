cmake_minimum_required (VERSION 2.6)
project (bamtools_test)

include_directories ("/usr/local/bamtools/include")

add_executable(bamtools_test test.cpp)