#!/usr/bin/env python

from optparse import OptionParser, OptionGroup
import re
import tempfile
from bs_align import output
from bs_align.bs_pair_end import *
from bs_align.bs_single_end import *
from bs_align.bs_rrbs import *
#import re
#from bs_utils.utils import *


if __name__ == '__main__':

    parser = OptionParser(usage="Usage: %prog {-i <single> | -1 <mate1> -2 <mate2>} -g <genome.fa> [options]")
    # option group 1
    opt_group = OptionGroup(parser, "For single end reads")
    opt_group.add_option("-i", "--input", type="string", dest="infilename",help="Input read file (FORMAT: sequences, qseq, fasta, fastq). Ex: read.fa or read.fa.gz", metavar="INFILE")
    parser.add_option_group(opt_group)

    # option group 2
    opt_group = OptionGroup(parser, "For pair end reads")
    opt_group.add_option("-1", "--input_1", type="string", dest="infilename_1",help="Input read file, mate 1 (FORMAT: sequences, qseq, fasta, fastq)", metavar="FILE")
    opt_group.add_option("-2", "--input_2", type="string", dest="infilename_2",help="Input read file, mate 2 (FORMAT: sequences, qseq, fasta, fastq)", metavar="FILE")
    opt_group.add_option("-I", "--minins",type = "int",dest = "min_insert_size", help="The minimum insert size for valid paired-end alignments [Default: %default]", default = 0)
    opt_group.add_option("-X", "--maxins",type = "int",dest = "max_insert_size", help="The maximum insert size for valid paired-end alignments [Default: %default]", default = 500)
    parser.add_option_group(opt_group)

    # option group 3
    opt_group = OptionGroup(parser, "Reduced Representation Bisulfite Sequencing Options")
    opt_group.add_option("-r", "--rrbs", action="store_true", dest="rrbs", default = False, help = 'Map reads to the Reduced Representation genome')
    opt_group.add_option("-c", "--cut-site", type="string",dest="cut_format", help="Cutting sites of restriction enzyme. Ex: MspI(C-CGG), Mael:(C-TAG), double-enzyme MspI&Mael:(C-CGG,C-TAG). [Default: %default]", metavar="pattern", default = "C-CGG")
    opt_group.add_option("-L", "--low", type = "int", dest="rrbs_low_bound", help="Lower bound of fragment length (excluding C-CGG ends) [Default: %default]", default = 20)
    opt_group.add_option("-U", "--up", type = "int", dest="rrbs_up_bound", help="Upper bound of fragment length (excluding C-CGG ends) [Default: %default]", default = 500)
    parser.add_option_group(opt_group)

    # option group 4
    opt_group = OptionGroup(parser, "General options")
    opt_group.add_option("-t", "--tag", type="string", dest="taginfo",help="[Y]es for undirectional lib, [N]o for directional [Default: %default]", metavar="TAG", default = 'N')
    opt_group.add_option("-s","--start_base",type = "int",dest = "cutnumber1", help="The first cycle of the read to be mapped [Default: %default]", default = 1)
    opt_group.add_option("-e","--end_base",type = "int",dest = "cutnumber2", help="The last cycle of the read to be mapped [Default: %default]", default = 200)
    opt_group.add_option("-a", "--adapter", type="string", dest="adapter_file",help="Input text file of your adaptor sequences (to be trimmed from the 3'end of the reads, ). "
                                                                                    "Input one seq for dir. lib., twon seqs for undir. lib. One line per sequence. "
                                                                                    "Only the first 10bp will be used", metavar="FILE", default = '')
    opt_group.add_option("--am",type = "int",dest = "adapter_mismatch", help="Number of mismatches allowed in adapter [Default: %default]", default = 0)
    opt_group.add_option("-g", "--genome", type="string", dest="genome",help="Name of the reference genome (should be the same as \"-f\" in bs_seeker2-build.py ) [ex. chr21_hg18.fa]")
    opt_group.add_option("-m", "--mismatches",type = "float", dest="no_mismatches",help="Number(>=1)/Percentage([0, 1)) of mismatches in one read. Ex: 4 (allow 4 mismatches) or 0.04 (allow 4% mismatches) [Default: %default]", default = 4)
    opt_group.add_option("--aligner", dest="aligner",help="Aligner program for short reads mapping: " + ', '.join(supported_aligners) + " [Default: %default]", metavar="ALIGNER", default = BOWTIE)
    opt_group.add_option("-p", "--path", dest="aligner_path", help="Path to the aligner program. Detected: " +' '*70+ '\t'.join(('%s: %s '+' '*70) % (al, aligner_path[al]) for al in sorted(supported_aligners)),
        metavar="PATH"
    )
    opt_group.add_option("-d", "--db", type="string", dest="dbpath",help="Path to the reference genome library (generated in preprocessing genome) [Default: %default]" , metavar="DBPATH", default = reference_genome_path)
    opt_group.add_option("-l", "--split_line",type = "int", dest="no_split",help="Number of lines per split (the read file will be split into small files for mapping. The result will be merged. [Default: %default]", default = 4000000)
    opt_group.add_option("-o", "--output", type="string", dest="outfilename",help="The name of output file [INFILE.bs(se|pe|rrbs)]", metavar="OUTFILE")
    opt_group.add_option("-f", "--output-format", type="string", dest="output_format",help="Output format: "+', '.join(output.formats)+" [Default: %default]", metavar="FORMAT", default = output.BAM)
    opt_group.add_option("--no-header", action="store_true", dest="no_SAM_header",help="Suppress SAM header lines [Default: %default]", default = False)
    opt_group.add_option("--temp_dir", type="string", dest="temp_dir",help="The path to your temporary directory [Detected: %default]", metavar="PATH", default = tempfile.gettempdir())
    opt_group.add_option("--XS",type = "string", dest="XS_filter",help="Filter definition for tag XS, format X,Y. X=0.8 and y=5 indicate that for one read, if #(mCH sites)/#(all CH sites)>0.8 and #(mCH sites)>5, then tag XS=1; or else tag XS=0. [Default: %default]", default = "0.5,5") # added by weilong

    opt_group.add_option("-M", "--multiple-hit", metavar="FileName", type="string", dest="Output_multiple_hit", default = None, help = 'File to store reads with multiple-hits')
    opt_group.add_option("-u", "--unmapped", metavar="FileName", type="string", dest="Output_unmapped_hit", default = None, help = 'File to store unmapped reads')

    opt_group.add_option("-v", "--version", action="store_true", dest="version",help="show version of BS-Seeker2", metavar="version", default = False)

    parser.add_option_group(opt_group)

    # option group 5
    opt_group = OptionGroup(parser, "Aligner Options",
        "You may specify any additional options for the aligner. You just have to prefix them with " +
        ', '.join('%s for %s' % (aligner_options_prefixes[aligner], aligner) for aligner in supported_aligners)+
        ', and BS Seeker will pass them on. For example: --bt-p 4 will increase the number of threads for bowtie to 4, '
        '--bt--tryhard will instruct bowtie to try as hard as possible to find valid alignments when they exist, and so on. '
        'Be sure that you know what you are doing when using these options! Also, we don\'t do any validation on the values.')
    parser.add_option_group(opt_group)


    #----------------------------------------------------------------
    # separate aligner options from BS Seeker options
    aligner_options = {}
    bs_seeker_options = []
    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        m = re.match(r'^%s' % '|'.join('(%s)'% aligner_options_prefixes[al] for al in supported_aligners), arg)
        if m:
            a_opt = arg.replace(m.group(0),'-',1)
            aligner_options[a_opt] = []
            while i + 1 < len(sys.argv) and sys.argv[i+1][0] != '-':
                aligner_options[a_opt].append(sys.argv[i+1])
                i += 1
            if len(aligner_options[a_opt]) == 0: # if it is a key-only option
                aligner_options[a_opt] = True
        else:
            bs_seeker_options.append(arg)
        i += 1

    (options, args) = parser.parse_args(args = bs_seeker_options)

    # if no options were given by the user, print help and exit
    if len(sys.argv) == 1:
        parser.print_help()
        exit(0)

    if options.version :
        show_version()
        exit (-1)
    else :
        show_version()

    # check parameters
    # input read files
    if options.infilename and (options.infilename_1 or options.infilename_2):
        error('-i and [-1|-2] options are exclusive. You should use only one of them.')

    if not (options.infilename or (options.infilename_1 and options.infilename_2)):
        error('You should set either -i or -1 and -2 options.')

    # Calculate the length of read
    if options.infilename :
        read_file = options.infilename
    elif options.infilename_1 :
        read_file = options.infilename_1
    else :
        error('You should at least specify -i or -1 options.')

    try :
        if read_file.endswith(".gz") : # support input file ending with ".gz"
            read_inf = gzip.open(read_file, "rb")
        else :
            read_inf=open(read_file,"r")
    except IOError :
        print "[Error] Cannot open input file : %s" % read_file
        exit(-1)
    oneline = read_inf.readline()
    oneline = read_inf.readline() # get the second line
    read_len = min(len(oneline), (options.cutnumber2-options.cutnumber1))
    read_inf.close()
    # mismatch allowed: bowtie 1,build-in parameter '-m'; bowtie 2, post-filter paramter
    # mismatch should no greater than the read length
    no_mismatches = float(options.no_mismatches)
    if (no_mismatches < 1) :
        int_no_mismatches=int(no_mismatches * read_len)
    else :
        int_no_mismatches=int(no_mismatches)

    str_no_mismatches=str(options.no_mismatches) # pass to specific mode


    # -t, directional / un-directional library
    asktag=str(options.taginfo).upper()
    if asktag not in 'YN':
        error('-t option should be either Y or N, not %s' % asktag)
    # -a
    if options.aligner not in supported_aligners:
        error('-a option should be: %s' % ' ,'.join(supported_aligners)+'.')
    # path for aligner
    aligner_exec = os.path.expanduser( os.path.join(options.aligner_path or aligner_path[options.aligner], options.aligner) )


    # -g
    if options.genome is None:
        error('-g is a required option')
    genome = os.path.split(options.genome)[1]
    genome_subdir = genome

    # try to guess the location of the reference genome for RRBS
    if options.rrbs:
        if options.rrbs_low_bound and options.rrbs_up_bound:
            if options.cut_format == "C-CGG" :
                genome_subdir += '_rrbs_%d_%d'  % (options.rrbs_low_bound, options.rrbs_up_bound)
            else :
                genome_subdir += '_rrbs_%s_%d_%d'  % ( re.sub(",","-",re.sub("-", "", options.cut_format)), options.rrbs_low_bound, options.rrbs_up_bound)
        else:
            possible_refs = filter(lambda dir: dir.startswith(genome+'_rrbs_'), os.listdir(options.dbpath))
            if len(possible_refs) == 1:
                genome_subdir = possible_refs[0]
            else:
                error('Cannot localize unambiguously the reference genome for RRBS. '
                      'Please, specify the options \"--low\" and \"--up\" that you used at the index-building step.\n'
                      'Possible choices are:\n' + '\n'.join([pr.split('_rrbs_')[-1].replace('_',', ') for pr in possible_refs]))

    db_path = os.path.expanduser(os.path.join(options.dbpath, genome_subdir + '_' + options.aligner))


    if not os.path.isdir(db_path):
        error('Index DIR \"' + genome_subdir + '..\" cannot be found in ' + options.dbpath +'.\n\tPlease run the bs_seeker2-build.py '
                            'to create it with the correct parameters for -g, -r, --low, --up and --aligner.')

    # default aligner options
    aligner_options_defaults = {
                                BOWTIE  : { '-e'              : 40*int_no_mismatches,
                                            '--nomaqround'    : True,
                                            '--norc'          : True,
                                            #'-k'              : 2,
                                            # -k=2; report two best hits, and filter by error rates
                                            '--quiet'         : True,
                                            '--best'          : True,
#                                            '--suppress'      : '2,5,6',
                                            '--sam'           : True,
                                            '--sam-nohead'    : True,
                                            '-p'              : 2
                                },
                                BOWTIE2 : {
                                            #'-M'              : 5,
                                            '--norc'          : True,
                                            '--quiet'         : True,
                                            '-p'              : 2,
                                            '--sam-nohead'    : True,
                                            # run bowtie2 in local mode by default
                                            '--local' : '--end-to-end' not in aligner_options,
                                            #'--mm'            : True,
                                            #'-k'              : 2
                                },
                                SOAP    : { '-v' : int_no_mismatches,
                                            '-p' : 2,
                                            '-r' : 2,
                                            '-M' : 4
                                          },
                                RMAP    : { '-M' : 2
                                            # to do # control for only mapping on + strand
                                          }

                                }

    if '--end-to-end' not in aligner_options:
        aligner_options_defaults[BOWTIE2].update({'-D' : 50})
        #aligner_options_defaults[BOWTIE2].update({'-D' : 50, '-R': 3, '-N': 0, '-L': 15, '-i' : 'S,1,0.50'})
    else:
        aligner_options_defaults[BOWTIE2].update({'-D' : 50, '-L': 15, '--score-min': 'L,-0.6,-0.6' })

    aligner_options = dict(aligner_options_defaults[options.aligner], **aligner_options)

    aligner_options_string = lambda : ' %s ' % (' '.join(opt_key +
                                                         (' ' + ' '.join(map(str,opt_val)) # join all values if the value is an array
                                                          if type(opt_val) is list else
                                                                ('' if type(opt_val) is bool and opt_val # output an empty string if it is a key-only option
                                                                 else ' ' +str(opt_val)) # output the value if it is a single value
                                                         )
                                                        for opt_key, opt_val in aligner_options.iteritems() if opt_val not in [None, False]))


#    tmp_path = (options.outfilename or options.infilename or options.infilename_1) +'-'+ options.aligner+ '-TMP'
#    clear_dir(tmp_path)

    options.output_format = options.output_format.lower()
    if options.output_format not in output.formats:
        error('Output format should be one of: ' + ', '.join(output.formats))

    if options.outfilename:
        outfilename = options.outfilename
        logfilename = outfilename
    elif options.infilename is not None:
        logfilename = options.infilename+'_'+ ('rr' if options.rrbs else '') + 'bsse'
        outfilename = logfilename + '.' + options.output_format
    else:
        logfilename = options.infilename_1+'_'+ ('rr' if options.rrbs else '') + 'bspe'
        outfilename = logfilename + '.' + options.output_format

    outfilename = os.path.expanduser(outfilename)
    logfilename = os.path.expanduser(logfilename)
    outfile = output.outfile(outfilename, options.output_format, deserialize(os.path.join(db_path, 'refname')), ' '.join(sys.argv), options.no_SAM_header)

    open_log(logfilename+'.bs_seeker2_log')

    aligner_title = options.aligner
    if options.aligner == BOWTIE2 :
        if '--end-to-end' in aligner_options :
            aligner_title = aligner_title + "-e2e"
        else:
            aligner_title = aligner_title + "-local"

    if options.aligner == BOWTIE :
        logm("Mode: Bowtie")
    elif options.aligner == BOWTIE2 :
        if '--end-to-end' not in aligner_options :
            logm("Mode: Bowtie2, local alignment")
        else :
            logm("Mode: Bowtie2, end-to-end alignment")


    tmp_path = tempfile.mkdtemp(prefix='bs_seeker2_%s_-%s-TMP-' % (os.path.split(outfilename)[1], aligner_title ), dir = options.temp_dir)


    (XS_x, XS_y) = options.XS_filter.split(",")
    XS_pct = float(XS_x)
    XS_count = int(XS_y)
    logm('Filter for tag XS: #(mCH)/#(all CH)>%.2f%% and #(mCH)>%d' % (XS_pct*100, XS_count))


    logm('Temporary directory: %s' % tmp_path)
    logm('Reduced Representation Bisulfite Sequencing: %s' % str(options.rrbs))
    if options.infilename is not None:
        logm('Single end')

        aligner_command = aligner_exec  + aligner_options_string() + \
                              { BOWTIE   : ' -k 2 %(reference_genome)s  -f %(input_file)s %(output_file)s',
                                BOWTIE2  : ' -k 2 -x %(reference_genome)s -f -U %(input_file)s -S %(output_file)s',
                                SOAP     : ' -D %(reference_genome)s.fa.index -o %(output_file)s -a %(input_file)s',
                                RMAP     : ' -c %(reference_genome)s.fa -o %(output_file)s %(input_file)s'
                              }[options.aligner]
        logm ('Aligner command: %s' % aligner_command)
        # single end reads
        if options.rrbs: # RRBS scan
            bs_rrbs(options.infilename,
                    asktag,
                    options.adapter_file,
                    options.cutnumber1,
                    options.cutnumber2,
                    options.no_split,
                    str_no_mismatches,
                    aligner_command,
                    db_path,
                    tmp_path,
                    outfile,
                    XS_pct,
                    XS_count,
                    options.adapter_mismatch,
                    options.Output_multiple_hit,
                    options.Output_unmapped_hit,
                    options.cut_format
                    )
        else: # Normal single end scan
            bs_single_end(  options.infilename,
                            asktag,
                            options.adapter_file,
                            options.cutnumber1,
                            options.cutnumber2,
                            options.no_split,
                            str_no_mismatches,
                            aligner_command,
                            db_path,
                            tmp_path,
                            outfile,
                            XS_pct,
                            XS_count,
                            options.adapter_mismatch,
                            options.Output_multiple_hit,
                            options.Output_unmapped_hit
                            )
    else:
        logm('Pair end')
        # pair end specific default options
        aligner_options = dict({BOWTIE: {'--fr'  : True,
                                         '-X'    : options.max_insert_size,
                                         '-I'    : options.min_insert_size if options.min_insert_size > 0 else None,
                                         '-a'    : True # "-k 2" in bowtie would not report the best two
                                },
                                BOWTIE2 : {
                                         '--fr'  : True,
                                         '-X'    : options.max_insert_size,
                                         '-I'    : options.min_insert_size if options.min_insert_size > 0 else None,
                                         '--no-discordant'  : True,
                                         '--no-mixed'       : True,
                                         '-k'    : 2
                                },
                                SOAP: {
                                        '-x' : options.max_insert_size,
                                        '-m' : options.min_insert_size if options.min_insert_size > 0 else 100
                                }}[options.aligner],
                               # integrating 'rmappe' is different from others
                                **aligner_options)

        aligner_command = aligner_exec + aligner_options_string() + \
                              { BOWTIE   : ' %(reference_genome)s  -f -1 %(input_file_1)s -2 %(input_file_2)s %(output_file)s',
                                BOWTIE2  : ' -x %(reference_genome)s  -f -1 %(input_file_1)s -2 %(input_file_2)s -S %(output_file)s',
                                SOAP     : ' -D %(reference_genome)s.fa.index -o %(output_file)s -a %(input_file_1)s -b %(input_file_2)s -2 %(output_file)s.unpaired' #,
                              #  RMAP     : #  rmappe, also paste two inputs into one file.
                             }[options.aligner]

        logm('Aligner command: %s' % aligner_command)

        if '--end-to-end' not in aligner_options:
            aligner_options_defaults[BOWTIE2].update({'-D' : 50})
        else:
            aligner_options_defaults[BOWTIE2].update({'-D' : 50, '-L': 15, '--score-min': 'L,-0.6,-0.6' })

        bs_pair_end(options.infilename_1,
                    options.infilename_2,
                    asktag,
                    options.adapter_file,
                    options.cutnumber1,
                    options.cutnumber2,
                    options.no_split,
                    str_no_mismatches,
                    aligner_command,
                    db_path,
                    tmp_path,
                    outfile,
                    XS_pct,
                    XS_count,
                    options.adapter_mismatch,
                    options.Output_multiple_hit,
                    options.Output_unmapped_hit
             )

    outfile.close()

