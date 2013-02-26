#!/usr/bin/python

import os
from optparse import OptionParser, OptionGroup
from bs_index.wg_build import *
from bs_index.rrbs_build import *
from bs_utils.utils import *


if __name__ == '__main__':

    parser = OptionParser()

    parser.add_option("-f", "--file", dest="filename",help="Input your reference genome file (fasta)", metavar="FILE")
    parser.add_option("--aligner", dest="aligner",help="Aligner program to perform the analysis: " + ', '.join(supported_aligners) + " [%default]", metavar="ALIGNER", default = BOWTIE2)
    parser.add_option("-p", "--path",   dest="aligner_path",help="Path to the aligner program. Defaults: " +' '*70+ '\t'.join(('%s: %s '+' '*70) % (al, aligner_path[al]) for al in sorted(supported_aligners)),
                  metavar="PATH")
    parser.add_option("-d", "--db", type="string", dest="dbpath",help="Path to the reference genome library (generated in preprocessing genome) [%default]", metavar="DBPATH", default = reference_genome_path)


    # RRBS options
    rrbs_opts = OptionGroup(parser, "Reduced Representation Bisulfite Sequencing Options",
                                "Use this options with conjuction of -r [--rrbs]")
    rrbs_opts.add_option("-r", "--rrbs", action="store_true", dest="rrbs", default = False, help = 'Preprocess the genome for analysis of Reduced Representation Bisulfite Sequencing experiments')
    rrbs_opts.add_option("-l", "--low",type= "int", dest="low_bound",help="lower bound [%default]", default = 50)
    rrbs_opts.add_option("-u", "--up", type= "int",dest="up_bound",help="upper bound [%default]", default = 300)
    parser.add_option_group(rrbs_opts)


    (options, args) = parser.parse_args()

    # if no options were given by the user, print help and exit
    if len(sys.argv) == 1:
        print parser.print_help()
        exit(0)

    rrbs = options.rrbs

    fasta_file=os.path.expanduser(options.filename)
    if fasta_file is None:
        error('Fasta file for the reference genome must be supported')

    if not os.path.isfile(fasta_file):
        error('%s cannot be found' % fasta_file)

    if options.aligner not in supported_aligners:
        error('-a option should be: %s' % ' ,'.join(supported_aligners)+'.')

    builder_exec = os.path.join(options.aligner_path or aligner_path[options.aligner],
                                {BOWTIE   : 'bowtie-build',
                                 BOWTIE2  : 'bowtie2-build',
                                 SOAP     : '2bwt-builder'}[options.aligner])

    build_command = builder_exec + { BOWTIE   : ' -f %(fname)s.fa %(fname)s',
                                     BOWTIE2  : ' -f %(fname)s.fa %(fname)s',
                                     SOAP     : ' %(fname)s.fa'
                                   }[options.aligner]


    print "Reference genome file: %s" % fasta_file
    print "Reduced Representation Bisulfite Sequencing: %s" % rrbs
    print "Builder path: %s" % builder_exec
    #---------------------------------------------------------------

    ref_path = options.dbpath

    if os.path.exists(ref_path):
        if not os.path.isdir(ref_path):
            error("%s must be a directory. Please, delete it or change the -d option." % ref_path)
    else:
        os.mkdir(ref_path)

    if rrbs: # RRBS preprocessing
        rrbs_build(fasta_file, build_command, ref_path, options.low_bound, options.up_bound, options.aligner)
    else: # Whole genome preprocessing
        wg_build(fasta_file, build_command, ref_path, options.aligner)
