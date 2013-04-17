#!/usr/bin/python

from optparse import OptionParser, OptionGroup
from bs_utils.utils import *

try :
    import pysam
except ImportError :
    print "[Error] Cannot import \"pysam\" package. Have you installed it?"
    exit(-1)

import gzip

def context_calling(seq, position):

    word=seq[position]
    word=word.upper()

    context="--"
    context_CH="--"
    if position + 2 < len(seq) and position - 2 >= 0:

        if word == "C":
            word2 = seq[position+1]
            context_CH = word + word2
            if word2 == "G":
                context = "CG"
            elif word2 in ['A','C','T']:
                word3 = seq[position+2]
                if word3 == "G":
                    context = "CHG"
                elif word3 in ['A','C','T']:
                    context="CHH"

        elif word == "G":
            word2 = seq[position-1]
            context_CH = word + word2
            context_CH = context_CH.translate(string.maketrans("ATCG", "TAGC"))
            if word2 == "C":
                context = "CG"
            elif word2 in ['A','G','T']:
                word3 = seq[position-2]
                if word3 == "C":
                    context = "CHG"
                elif word3 in ['A','G','T']:
                    context = "CHH"

    return word, context, context_CH



if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("-i", "--input", type="string", dest="infilename",help="BAM output from bs_seeker2-align.py", metavar="INFILE")
    parser.add_option("-d", "--db", type="string", dest="dbpath",help="Path to the reference genome library (generated in preprocessing genome) [Default: %default]" , metavar="DBPATH", default = reference_genome_path)
    parser.add_option("-o", "--output-prefix", type="string", dest="output_prefix",help="The output prefix to create ATCGmap and wiggle files [INFILE]", metavar="OUTFILE")

    parser.add_option("--wig", type="string", dest="wig_file",help="The output .wig file [INFILE.wig]", metavar="OUTFILE")
    parser.add_option("--CGmap", type="string", dest="CGmap_file",help="The output .CGmap file [INFILE.CGmap]", metavar="OUTFILE")
    parser.add_option("--ATCGmap", type="string", dest="ATCGmap_file",help="The output .ATCGmap file [INFILE.ATCGmap]", metavar="OUTFILE")

    parser.add_option("-x", "--rm-SX", action="store_true", dest="RM_SX",help="Removed reads with tag \'XS:i:1\', which would be considered as not fully converted by bisulfite treatment [Default: %default]", default = False)

    parser.add_option("-r", "--read-no",type = "int", dest="read_no",help="The least number of reads covering one site to be shown in wig file [Default: %default]", default = 1)
    parser.add_option("-v", "--version", action="store_true", dest="version",help="show version of BS-Seeker2", metavar="version", default = False)

    (options, args) = parser.parse_args()


    # if no options were given by the user, print help and exit
    if len(sys.argv) == 1:
        print parser.print_help()
        exit(0)

    if options.version :
        show_version()
        exit (-1)
    else :
        show_version()


    if options.infilename is None:
        error('-i option is required')
    if not os.path.isfile(options.infilename):
        error('Cannot find input file: %s' % options.infilename)

    open_log(options.infilename+'.call_methylation_log')
    db_d = lambda fname:  os.path.join( os.path.expanduser(options.dbpath), fname) # bug fixed, weilong

    logm('sorting BS-Seeker alignments')
    sorted_input_filename = options.infilename+'_sorted'
    pysam.sort(options.infilename, sorted_input_filename)
    sorted_input_filename += '.bam'
    logm('indexing sorted alignments')
    pysam.index(sorted_input_filename)

    logm('calculating methylation levels')
    ATCGmap_fname = options.ATCGmap_file or ((options.output_prefix or options.infilename) + '.ATCGmap.gz')
    ATCGmap = gzip.open(ATCGmap_fname, 'wb')

    CGmap_fname = options.CGmap_file or ((options.output_prefix or options.infilename) + '.CGmap.gz')
    CGmap = gzip.open(CGmap_fname, 'wb')

    wiggle_fname = options.wig_file or ((options.output_prefix or options.infilename) + '.wig')
    wiggle = open(wiggle_fname, 'w')

    sorted_input = pysam.Samfile(sorted_input_filename, 'rb')
    
    chrom = None
    nucs = ['A', 'T', 'C', 'G', 'N']
    ATCG_fwd = dict((n, 0) for n in nucs)
    ATCG_rev = dict((n, 0) for n in nucs)
    for col in sorted_input.pileup():
        col_chrom = sorted_input.getrname(col.tid)
        if chrom != col_chrom:
            chrom = col_chrom
            chrom_seq = deserialize(db_d(chrom))
            wiggle.write('variableStep chrom=%s\n' % chrom)

        for n in nucs:
            ATCG_fwd[n] = 0
            ATCG_rev[n] = 0

        nuc, context, subcontext = context_calling(chrom_seq, col.pos)
        total_reads = 0



        for pr in col.pileups:
        #     print pr
             if (not pr.indel) : # skip indels
                #if ( (options.RM_SX) and (pr.alignment.tags[1][1] == 1) ):
                ##=== Fixed error reported by Roberto
                #print options.RM_SX,  dict(pr.alignment.tags)["XS"]
                if ( (options.RM_SX) and (dict(pr.alignment.tags)["XS"] == 1) ):
                    # print "Debug: ", options.RM_SX, pr.alignment.tags[1]
                    # when need to filter and read with tag (XS==1), then remove the reads
                    continue

                if pr.qpos >= len(pr.alignment.seq):
                    print 'WARNING: read %s has an invalid alignment. Discarding.. ' % pr.alignment.qname
                    continue
                read_nuc = pr.alignment.seq[pr.qpos]
         #       print "read_nuc=", read_nuc
                if pr.alignment.is_reverse:
                    ATCG_rev[read_nuc] += 1
                else:
                    ATCG_fwd[read_nuc] += 1

                if read_nuc != 'N':
                    total_reads += 1

        cnts = lambda d: '\t'.join(str(d[n]) for n in nucs)
        fwd_counts = cnts(ATCG_fwd)
        rev_counts = cnts(ATCG_rev)

        meth_level = None
        meth_cytosines = 0
        unmeth_cytosines = 0

        if nuc == 'C':
            # plus strand: take the ratio of C's to T's from reads that come from the forward strand
            meth_cytosines = ATCG_fwd['C']
            unmeth_cytosines = ATCG_fwd['T']

        elif nuc == 'G':
            # minus strand: take the ratio of G's to A's from reads that come from the reverse strand
            meth_cytosines = ATCG_rev['G']
            unmeth_cytosines = ATCG_rev['A']

        if meth_cytosines + unmeth_cytosines > 0:
            meth_level = float(meth_cytosines)/(meth_cytosines + unmeth_cytosines)

        pos = col.pos + 1

        meth_level_string = str(meth_level) if meth_level is not None else 'na'
        ATCGmap.write('%(chrom)s\t%(nuc)s\t%(pos)d\t%(context)s\t%(subcontext)s\t%(fwd_counts)s\t%(rev_counts)s\t%(meth_level_string)s\n' % locals())
#
        all_cytosines = meth_cytosines + unmeth_cytosines 
        if (meth_level is not None) and (all_cytosines >= options.read_no):
        #    print all_cytosines
            if nuc == 'C':
                wiggle.write('%d\t%f\n' % (pos, meth_level))
            else :
                wiggle.write('%d\t-%f\n' % (pos, meth_level))
            CGmap.write('%(chrom)s\t%(nuc)s\t%(pos)d\t%(context)s\t%(subcontext)s\t%(meth_level_string)s\t%(meth_cytosines)s\t%(all_cytosines)s\n' % locals())

    logm('Wiggle: %s'% wiggle_fname)
    logm('ATCGMap: %s' % ATCGmap_fname)
    logm('CGmap: %s' % CGmap_fname)

