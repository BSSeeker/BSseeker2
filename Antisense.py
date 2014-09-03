#!/usr/bin/env python
# 2014-08-13

import sys
#import os
#import os.path
#import re
import gzip



# ===========================================
bc_dict = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N',
           'a':'t', 'c':'g', 'g':'c', 't':'a', 'n':'n'}

def bc ( ch ) :
    if ch in bc_dict :
        return bc_dict[ch]
    else :
        return ch

def AntisenseRead ( read ) :
    read="".join([ bc(i) for i in read])
    return read[::-1]

# 
def AntisenseFastq (infile, line1) :
    try:
        if infile :
            if infile.endswith(".gz") :
                IN = gzip.open(infile, 'rb')
            else :
                IN = open(infile, 'r')
        else :
            IN = sys.stdin
    except IOError:
        print("\n[Error]:\n\t File cannot be open: %s" % infile)
        exit(-1)


    count = 0
    total = 0
    first = True
    while 1 :
        if infile or (not first) :
            line1 = IN.readline().strip()
        else :
            first = False
        if not line1 :
            break

        line2 = IN.readline().strip()
        line3 = IN.readline().strip()
        line4 = IN.readline().strip()

        print(line1)
        print(AntisenseRead(line2))
        print(line3)
        print(line4[::-1])

    if infile :
        IN.close()



# ===============================
#

def AntisenseSequence (infile, line) :
    try:
        if infile :
            if infile.endswith(".gz") :
                IN = gzip.open(infile, 'rb')
            else :
                IN = open(infile, 'r')
        else :
            IN = sys.stdin
    except IOError:
        print("\n[Error]:\n\t File cannot be open: ", infile)
        exit(-1)

    if infile is None :
        print( AntisenseRead(line) )

    for line in IN :
        line = line.strip()
        print( AntisenseRead(line) )

    if infile :
        IN.close()


# ===============================
# SN603   WA047   6       1101    41.40   99.10   0       1       .GGGA.......TTAG..............       @SZ__@@@@@@@RR__@@@@@@@@@@@@@@       0
#
def AntisenseQseq (infile, line, keepquality) :
    try:
        if infile :
            if infile.endswith(".gz") :
                IN = gzip.open(infile, 'rb')
            else :
                IN = open(infile, 'r')
        else :
            IN = sys.stdin
    except IOError:
        print("\n[Error]:\n\t File cannot be open: %s" % infile)
        exit(-1)

    if infile is None :
        tokens = line.split()
        if ( (not keepquality and tokens[10] == "1") or keepquality ) :
            tokens[8]=AntisenseRead(tokens[8])
            tokens[9]=tokens[9][::-1]
            print( "\t".join(tokens) )

    for line in IN :
        line = line.strip()
        tokens = line.split()
        if ( (not keepquality and tokens[10] == "1") or keepquality ) :
            tokens[8]=AntisenseRead(tokens[8])
            tokens[9]=tokens[9][::-1]
            print( "\t".join(tokens) )

    if infile :
        IN.close()

# ===============================
#
def printFasta( read, length ) :
    T_len = len(read)
    length=int(length)
    #print "printFasta"
    if length >0 :
        for start in xrange(0, T_len, length) :
            if start+length > T_len :
                end = T_len
            else :
                end = start+length
            print(read[start:end])
    else :
        print(read)


def AntisenseFasta (infile, line, length) :
    try:
        if infile :
            if infile.endswith(".gz") :
                IN = gzip.open(infile, 'rb')
            else :
                IN = open(infile, 'r')
        else :
            IN = sys.stdin

    except IOError:
        print("\n[Error]:\n\t File cannot be open: ", infile)
        exit(-1)

    name = ""
    read = ""

    if infile is None :
        read = ""

    for line in IN :
        line = line.strip()
        if (line[0] == '>') :
            if ( name != "" ) :
                printFasta( AntisenseRead(read), length )
            print(line)
            name=line
            read=""
        else :
            read = read + line
    if (name != "") :
            printFasta( AntisenseRead(read), length )




# ==============================
#  Decide branch functions for different format

def Antisense (infile, keepquality, length):
    try:
        if infile :
            if infile.endswith(".gz") :
                IN = gzip.open(infile, 'rb')
            else :
                IN = open(infile, 'r')
        else :
            IN = sys.stdin

    except IOError:
        print("\n[Error]:\n\t File cannot be open: ", infile)
        exit(-1)
    i = 0

    line = IN.readline().strip()
    tokens = line.split()
    input_format = ""
    if line[0]=="@" : 
        input_format = "FastQ"
    elif len(tokens) == 1 and line[0] != ">":
        input_format = "sequence"
    elif len(tokens) == 11:
        input_format = "qseq"
    elif line[0]==">" :
        input_format = "fasta"

    if infile :
        IN.close()

    if input_format == "FastQ" :
        AntisenseFastq(infile, line)
    elif input_format == "sequence" :
        AntisenseSequence(infile, line)
    elif input_format == "qseq" :
        AntisenseQseq(infile, line, keepquality)
    elif input_format == "fasta" :
        AntisenseFasta(infile, line, length)


from optparse import OptionParser

# ===========================================
def main():
    usage = "Usage: %prog [-i <in>] [-o <output>] [-l 50]\n" \
            "Author : Guo, Weilong; guoweilong@gmail.com; 2012-11-10\n" \
            "Last Update: 2014-08-13\n" \
            "Description: Unique reads for qseq/fastq/fasta/sequencce,\n" \
            "    Discription: To get the revesed complementary sequences of the input.fa.\n" \
	    "    The upper/lower cases are kept for corresponding sites.\n" \
            "Example:\n" \
            "input.fa:  ACCGTTCCTTG\n" \
            "output.fa: CAAGGAACGGT"
    parser = OptionParser(usage)
    parser.add_option("-i", dest="infile",
                  help="input file, format: qseq/fastq/fasta/sequence", metavar="FILE")
    parser.add_option("-o", dest="outfile",
                  help="Name of the output file", metavar="FILE")
    parser.add_option("-l", dest="length", default=0,
                  help="Length to show in one line for fasta. [Default: show all in one line]", metavar="INT")
    parser.add_option("-k", dest="keepquality", default = False, action = "store_true",
                  help="Would not filter low quality reads if specified")
    (options, args) = parser.parse_args()
    
    if (options.outfile is not None) :
        if options.outfile.endswith('.gz') :
            sys.stdout = gzip.open(options.outfile, 'wb')
        else :
            sys.stdout = open(options.outfile, 'w')

    Antisense (options.infile, options.keepquality, options.length)

# ===========================================
if __name__ == "__main__":
    main()

