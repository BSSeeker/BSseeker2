#!/usr/bin/env python

import os
from optparse import OptionParser, OptionGroup
from bs_index.wg_build import *
from bs_index.rrbs_build import *
from bs_utils.utils import *


if __name__ == '__main__':

    parser = OptionParser()

    parser.add_option("-f", "--file", dest="filename", help="Input your reference genome file (fasta)", metavar="FILE")
    parser.add_option("--aligner", dest="aligner", help="Aligner program to perform the analysis: " + ', '.join(supported_aligners) + " [Default: %default]", metavar="ALIGNER", default = BOWTIE)
    parser.add_option("-p", "--path", dest="aligner_path", help="Path to the aligner program. Detected: " +' '*70+ '\t'.join(('%s: %s '+' '*70) % (al, aligner_path[al]) for al in sorted(supported_aligners)),
                  metavar="PATH")
    parser.add_option("-d", "--db", type="string", dest="dbpath", help="Path to the reference genome library (generated in preprocessing genome) [Default: %default]", metavar="DBPATH", default = reference_genome_path)

    parser.add_option("-v", "--version", action="store_true", dest="version", help="show version of BS-Seeker2", default=False)

    # RRBS options
    rrbs_opts = OptionGroup(parser, "Reduced Representation Bisulfite Sequencing Options",
                                "Use this options with conjuction of -r [--rrbs]")
    rrbs_opts.add_option("-r", "--rrbs", action="store_true", dest="rrbs", help = 'Build index specially for Reduced Representation Bisulfite Sequencing experiments. Genome other than certain fragments will be masked. [Default: %default]', default = False)
    rrbs_opts.add_option("-l", "--low",type= "int", dest="low_bound", help="lower bound of fragment length (excluding recognition sequence such as C-CGG) [Default: %default]", default = 20)
    rrbs_opts.add_option("-u", "--up", type= "int", dest="up_bound", help="upper bound of fragment length (excluding recognition sequence such as C-CGG ends) [Default: %default]", default = 500)
    rrbs_opts.add_option("-c", "--cut-site", type= "string", dest="cut_format", help="Cut sites of restriction enzyme. Ex: MspI(C-CGG), Mael:(C-TAG), double-enzyme MspI&Mael:(C-CGG,C-TAG). [Default: %default]", default = "C-CGG")
    parser.add_option_group(rrbs_opts)


    (options, args) = parser.parse_args()

    # if no options were given by the user, print help and exit
    if len(sys.argv) == 1:
        parser.print_help()
        exit(0)

    if options.version :
        show_version()
        exit (-1)
    else :
        show_version()

    rrbs = options.rrbs

    if options.filename is not None :
        fasta_file=os.path.expanduser(options.filename)
    else :
        error("Please specify the genome file (Fasta) using \"-f\"")

    if fasta_file is None:
        error('Fasta file for the reference genome must be supported')

    if not os.path.isfile(fasta_file):
        if os.path.isfile(os.path.join(os.dbpath, fasta_file)):
            # Search for os.dbpath to check if the genome file is stored there.
            fasta_file = os.path.join(os.dbpath, fasta_file)
        else:
            error('%s cannot be found' % fasta_file)

    if options.aligner not in supported_aligners:
        error('-a option should be: %s' % ' ,'.join(supported_aligners)+'.')

    builder_exec = os.path.join(options.aligner_path or aligner_path[options.aligner],
                                {BOWTIE   : 'bowtie-build',
                                 BOWTIE2  : 'bowtie2-build',
                                 SOAP     : '2bwt-builder',
                                 RMAP     : '' # do nothing
                                }[options.aligner])

    build_command = builder_exec + { BOWTIE   : ' -f %(fname)s.fa %(fname)s',
                                     BOWTIE2  : ' -f %(fname)s.fa %(fname)s',
                                     SOAP     : ' %(fname)s.fa'
                                   }[options.aligner]


    print "Reference genome file: %s" % fasta_file
    print "Reduced Representation Bisulfite Sequencing: %s" % rrbs
    print "Short reads aligner you are using: %s" % options.aligner
    print "Builder path: %s" % builder_exec

    #---------------------------------------------------------------

    if not os.path.isfile( builder_exec ) :
        error("Cannot file program %s for execution." % builder_exec)

    ref_path = options.dbpath

    if os.path.exists(ref_path):
        if not os.path.isdir(ref_path):
            error("%s must be a directory. Please, delete it or change the -d option." % ref_path)
    else:
        os.mkdir(ref_path)

    if rrbs: # RRBS preprocessing
        rrbs_build(fasta_file, build_command, ref_path, options.low_bound, options.up_bound, options.aligner, options.cut_format)
    else: # Whole genome preprocessing
        wg_build(fasta_file, build_command, ref_path, options.aligner)

