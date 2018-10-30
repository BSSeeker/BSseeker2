#!/usr/bin/env sh
# Guo, Weilong; guoweilong@126.com; 2017-10-26

# This file is only useful when the genome has lots of contigs
# That you want to threading these contigs as a few chromesomes for mapping
# The you need to Change the coodinate back to the original position.

g++ ChangeCoodinate.cpp -o ChangeCoodinate
g++ ThreadFasta.cpp  -o ThreadFasta


