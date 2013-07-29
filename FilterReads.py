#!/usr/bin/env python

import sys
import os
import os.path
import re

# Define the Dictionary & the Reverse Dictionary for barcode

# ===========================================
# 
def FilterFastq (infile, outfile) :
    try:
        INFILE = open(infile, 'r')
    except IOError:
        print "\n[Error]:\n\t File cannot be open: ", infile;
        exit(-1)
    try:
        OUTFILE = open(outfile, 'w')
    except IOError:
        print "\n[Error]:\n\t File cannot be open: ", outfile;
        exit(-1)
    Dict = {}
    count = 0
    total = 0
    while 1 :
        line1 = INFILE.readline()
        if not line1 :
            break
        line1 = line1
        line2 = INFILE.readline()
        line3 = INFILE.readline()
        line4 = INFILE.readline()
        total += 1
 #       print line2
        if ( line2 in Dict ) :
            Dict[line2] += 1
        else :
            Dict[line2] = 1
            OUTFILE.write(line1)
            OUTFILE.write(line2)
            OUTFILE.write(line3)
            OUTFILE.write(line4)
            count += 1
    INFILE.close()
    OUTFILE.close()
    print "Count of total raw reads: ", total
    print "Count of reads left: ", count
    print "Percentage of kept reads: %.2f" % (count * 100.0 /total), "%"


# ===============================
#

def FilterSequence (infile, outfile) :
    try:
        INFILE = open(infile, 'r')
    except IOError:
        print "\n[Error]:\n\t File cannot be open: ", infile;
        exit(-1)
    try:
        OUTFILE = open(outfile, 'w')
    except IOError:
        print "\n[Error]:\n\t File cannot be open: ", outfile;
        exit(-1)
    Dict = {}
    count = 0
    total = 0
    for line in INFILE:
        line = line.strip()
        total += 1
        if ( line in Dict ) :
            Dict[line] += 1
        else :
            Dict[line] = 1
            count += 1
            OUTFILE.write(line + "\n")
    INFILE.close()
    OUTFILE.close()
    print "Count of total raw reads: ", total
    print "Count of reads left: ", count
    print "Percentage of kept reads: %.2f" % (count * 100.0 /total), "%"

# ===============================
# SN603   WA047   6       1101    41.40   99.10   0       1       .GGGA.......TTAG..............       @SZ__@@@@@@@RR__@@@@@@@@@@@@@@       0
#
def FilterQseq (infile, outfile, keepquality) :
    if keepquality :
        print "User specified '-k' and read quality will not be considered"
    else :
        print "Reads with PF=0 will be filtered"

    try:
        INFILE = open(infile, 'r')
    except IOError:
        print "\n[Error]:\n\t File cannot be open: ", infile;
        exit(-1)
    try:
        OUTFILE = open(outfile, 'w')
    except IOError:
        print "\n[Error]:\n\t File cannot be open: ", outfile;
        exit(-1)
    Dict = {}
    count = 0
    total = 0
    for line in INFILE:
        tokens = line.strip().split()
        total += 1
        if ( (not keepquality and tokens[10] == "1") or keepquality ) :
            if ( tokens[8] in Dict ) :
                Dict[tokens[8]] += 1
            else :
                Dict[tokens[8]] = 1
                count += 1
                OUTFILE.write(line)

    INFILE.close()
    OUTFILE.close()
    print "Count of total raw reads: ", total
    print "Count of reads left: ", count
    print "Percentage of kept reads: %.2f" % (count * 100.0 /total), "%"


# ===============================
#

def FilterFasta (infile, outfile) :
    try:
        INFILE = open(infile, 'r')
    except IOError:
        print "\n[Error]:\n\t File cannot be open: ", infile;
        exit(-1)
    try:
        OUTFILE = open(outfile, 'w')
    except IOError:
        print "\n[Error]:\n\t File cannot be open: ", outfile;
        exit(-1)
    Dict = {}
    count = 0
    total = 0
    name = ""
    read = ""
    for line in INFILE:
        line = line.strip()
        if (line[0] == '>') :
            total += 1
            if ( name != "" ) :
                if ( read in Dict ) :
                    Dict[read] += 1
                else :
                    Dict[read] = 1
                    count += 1
                    OUTFILE.write(">" + name + "\n")
                    OUTFILE.write(read + "\n")
            name = line[1:]
            read = ""
        else :
            read = read + line
    if (name != "") :
        if ( read in Dict ) :
            Dict[read] += 1
        else :
            Dict[read] = 1
            count += 1
            OUTFILE.write(">" + name + "\n")
            OUTFILE.write(read + "\n")
    print "Count of total raw reads: ", total
    print "Count of reads left: ", count
    print "Percentage of kept reads: %.2f" % (count * 100.0 /total), "%"



# ==============================
#  Decide branch functions for different format

def FilterReads (infile, outfile, keepquality):
    try:
        INFILE = open(infile, 'r')
    except IOError:
        print "\n[Error]:\n\t File cannot be open: ", infile;
        exit(-1)
    i = 0
    # Count the numbers of barcodes in first 10000 lines
    line = INFILE.readline()
    tokens = line.strip().split()
    input_format="";
    if line[0]=="@" : 
        input_format = "FastQ"
        n_fastq = 0
    elif len(tokens) == 1 and line[0] != ">":
        input_format = "sequence"
    elif len(tokens) == 11:
        input_format = "qseq"
    elif line[0]==">" :
        input_format = "fasta"
    INFILE.close()

    print "Input file format detected: ", input_format

    if input_format == "FastQ" :
        FilterFastq(infile, outfile)
    elif input_format == "sequence" :
        FilterSequence(infile, outfile)
    elif input_format == "qseq" :
        FilterQseq(infile, outfile, keepquality)
    elif input_format == "fasta" :
        FilterFasta(infile, outfile)


from optparse import OptionParser

# ===========================================
def main():
    usage = "Usage: %prog -i <input> -o <output> [-k]\n" \
            "Author : Guo, Weilong; guoweilong@gmail.com; 2012-11-10\n" \
            "Last Update: 2013-04-01" \
            "Description: Unique reads for qseq/fastq/fasta/sequencce,\n" \
            "       and filter low quality reads in qseq file."
    parser = OptionParser(usage)
    parser.add_option("-i", dest="infile", 
                  help="Name of the input qseq/fastq/fasta/sequence file", metavar="FILE")
    parser.add_option("-o", dest="outfile",
                  help="Name of the output file", metavar="FILE")
    parser.add_option("-k", dest="keepquality", default = False, action = "store_true",
                  help="Would not filter low quality reads if specified")
    (options, args) = parser.parse_args()
    
    if (options.infile is None) or (options.outfile is None) :
        parser.print_help()
        exit(-1)
    FilterReads(options.infile, options.outfile, options.keepquality)


# ===========================================
if __name__ == "__main__":
    main()

