import gzip
import os

import datetime
import re
import shlex
import shutil
import string
import subprocess
import types
#from itertools import izip

import marshal
import sys



_rc_trans = string.maketrans('ACGT', 'TGCA')
def reverse_compl_seq(strseq):
    return strseq.translate(_rc_trans)[::-1]



def show_version() :
    print ""
    print "     BS-Seeker2 v2.1.5 - Dec. 21, 2017"
    print ""



"""
IUPAC nucleotide code Base
    A	Adenine
    C	Cytosine
    G	Guanine
    T (or U)	Thymine (or Uracil)
    R	A or G
    Y	C or T
    S	G or C
    W	A or T
    K	G or T
    M	A or C
    B	C or G or T
    D	A or G or T
    H	A or C or T
    V	A or C or G
    N	any base
    . or -	gap
"""


def IUPAC ( nuc ) :
    if nuc == 'R' :
        return ('A','G')
    elif nuc == 'Y' :
        return ('C', 'T')
    elif nuc == 'S' :
        return ('G', 'C')
    elif nuc == 'W' :
        return ('A', 'T')
    elif nuc == 'K' :
        return ('G','T')
    elif nuc == 'M' :
        return ('A','C')
    elif nuc == 'B' :
        return ('C', 'G', 'T')
    elif nuc == 'D' :
        return ('A', 'G', 'T')
    elif nuc == 'H' :
        return ('A', 'C', 'T')
    elif nuc == 'V' :
        return ('A', 'C', 'G')
    elif nuc == 'N' :
        return ('A', 'C', 'G', 'T')
    else :
        return (nuc)


def uniq(inlist):
    # order preserving
    uniques = []
    for item in inlist:
        if item not in uniques:
            uniques.append(item)
    return uniques

from itertools import product

def EnumerateIUPAC ( context_lst ) :
    tag_list = []
#    context_lst = [context]
    for one_context in context_lst :
        for m in product(*[ IUPAC(i) for i in list(one_context)]) :
            tag_list.append(''.join(m))
    return uniq(tag_list)

from itertools import product

# example: cut3_context="CGG"
# return generator for : ["CGG","TGG"]
# wild-card C to both C and T
def Enumerate_C_to_CT ( cut3_context_lst ) :
    tag_list = []
    for context in cut3_context_lst :
        for m in product(*[i if (i is not 'C') else ('C','T') for i in context]) :
            tag_list.append(''.join(m))
    return uniq(tag_list)

#-------------------------------------------------------------------------------------

# set a reasonable defaults
def find_location(program):
    def is_exe(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)
    for path in os.environ["PATH"].split(os.pathsep):
        if is_exe(os.path.join(path, program)):
            return path

    return None

BOWTIE = 'bowtie'
BOWTIE2 = 'bowtie2'
SOAP = 'soap'
RMAP = 'rmap'

supported_aligners = [
                      BOWTIE,
                      BOWTIE2,
                      SOAP,
                      RMAP
                    ]

aligner_options_prefixes = { BOWTIE  : '--bt-',
                             BOWTIE2 : '--bt2-',
                             SOAP    : '--soap-',
                             RMAP    : '--rmap-' }

#aligner_path = dict((aligner, os.path.expanduser(find_location(aligner) or default_path))
#                    for aligner, default_path in
#                           [(BOWTIE,'~/bowtie/'),
#                            (BOWTIE2, '~/bowtie2/'),
#                            (SOAP, '~/soap/'),
#                            (RMAP, '~/rmap/bin')
#                            ])
aligner_path = dict((aligner, os.path.expanduser(find_location(aligner) or "None"))
                    for aligner in [(BOWTIE), (BOWTIE2), (SOAP), (RMAP)])

reference_genome_path = os.path.join(os.path.split(globals()['__file__'])[0],'reference_genomes')



def error(msg):
    print >> sys.stderr, 'ERROR: %s' % msg
    exit(1)


global_stime = datetime.datetime.now()
def elapsed(msg = None):
    print "[%s]" % msg if msg is not None else "+", "Last:" , datetime.datetime.now() - elapsed.stime, '\tTotal:', datetime.datetime.now() - global_stime

    elapsed.stime = datetime.datetime.now()

elapsed.stime = datetime.datetime.now()


def open_log(fname):
    open_log.logfile = open(fname, 'w', 1)

def logm(message):
    log_message = "[%s] %s\n" % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), message)
    print log_message,
    open_log.logfile.write(log_message)

def close_log():
    open_log.logfile.close()




def clear_dir(path):
    """ If path does not exist, it creates a new directory.
        If path points to a directory, it deletes all of its content.
        If path points to a file, it raises an exception."""

    if os.path.exists(path):
        if not os.path.isdir(path):
            error("%s is a file. Please, delete it manually!" % path)
        else:
            for the_file in os.listdir(path):
                file_path = os.path.join(path, the_file)
                try:
                    if os.path.isfile(file_path):
                        os.unlink(file_path)
                    elif os.path.isdir(file_path):
                        shutil.rmtree(file_path)
                except Exception, e:
                    print e
    else:
        os.mkdir(path)


def delete_files(*filenames):
#    return
    """ Deletes a number of files. filenames can contain generator expressions and/or lists, too"""

    for fname in filenames:
        if type(fname) in [list, types.GeneratorType]:
            delete_files(*list(fname))
        elif os.path.isdir(fname):
            shutil.rmtree(fname)
        else:
            os.remove(fname)

def split_file(filename, output_prefix, nlines):
    """ Splits a file (equivalend to UNIX split -l ) """
    fno = 0
    lno = 0
    if filename.endswith(".gz") :
        INPUT = gzip.open(filename, 'rb')
    else :
        INPUT = open(filename, 'r')
    output = None
    for l in INPUT:
        if lno == 0:
            fno += 1
            if output is not None: output.close()
            output = open('%s%d' % (output_prefix, fno), 'w')
            lno = nlines
        output.write(l)
        lno -= 1
    output.close()
    INPUT.close()

def isplit_file(filename, output_prefix, nlines):
    """ Splits a file (equivalend to UNIX split -l ) """
    fno = 0
    lno = 0
    try :
        input = (gzip.open if filename.endswith('.gz') else open)(filename, 'r')
    except IOError :
        print "[Error] Cannot find file : %s !" % filename
        exit(-1)
    output = None
    output_fname = None
    for l in input:
        if lno == 0:

            if output is not None:
                output.close()
                yield output_fname

            fno += 1
            output_fname = '%s%d' % (output_prefix, fno)
#            output_fname = '%s_0' % output_prefix
            output = open(output_fname, 'w')
            lno = nlines
        output.write(l)
        lno -= 1
    if output is not None:
        output.close()
    yield output_fname

    input.close()

def read_fasta(fasta_file):
    """ Iterates over all sequences in a fasta file. One at a time, without reading the whole file into the main memory.
    """

    try :
        input = (gzip.open if fasta_file.endswith('.gz') else open)(fasta_file)
    except IOError:
        print "[Error] Cannot find fasta file : %s !" % fasta_file
        exit(-1)
    sanitize = re.compile(r'[^ACTGN]')
    sanitize_seq_id = re.compile(r'[^A-Za-z0-9]')

    chrom_seq = []
    chrom_id = None
    seen_ids = set()

    for line in input:
        if line[0] == '>':
            if chrom_id is not None:
                yield chrom_id, ''.join(chrom_seq)

            chrom_id = sanitize_seq_id.sub('_', line.split()[0][1:])

            if chrom_id in seen_ids:
                error('BS Seeker found identical sequence ids (id: %s) in the fasta file: %s. Please, make sure that all sequence ids are unique and contain only alphanumeric characters: A-Za-z0-9_' % (chrom_id, fasta_file))
            seen_ids.add(chrom_id)

            chrom_seq = []

        else:
            chrom_seq.append(sanitize.sub('N', line.strip().upper()))

    yield chrom_id, ''.join(chrom_seq)

    input.close()


def serialize(obj, filename):
    """ Be careful what you serialize and deseriazlize! marshal.load is not secure!
    """
    output = open(filename+'.data', 'w')
    marshal.dump(obj, output)
    output.close()

def deserialize(filename):
    """ Be careful what you serialize and deseriazlize! marshal.load is not secure!
    """
    try:
        input = open(filename+'.data')
    except IOError:
        print "\n[Error]:\n\t Cannot find file: %s.data" % filename
        exit(-1)

    obj = marshal.load(input)
    input.close()
    return obj   



def run_in_parallel(commands):

    commands = [(cmd[0], open(cmd[1], 'w')) if type(cmd) is tuple else (cmd, None) for cmd in commands]

    #logm('Starting commands:\n' + '\n'.join([cmd for cmd, stdout in commands]))
    logm('Starting commands:')
    for cmd, stdout in commands :
        logm("Launched: "+cmd)
    for i, proc in enumerate([subprocess.Popen(args = shlex.split(cmd), stdout = stdout) for cmd, stdout in commands]):
        return_code = proc.wait()
        logm('Finished: ' + commands[i][0])
        if return_code != 0:
            error('%s \nexit with an error code: %d. Please, check the log files.' % (commands[i], return_code))
    for _, stdout in commands:
        if stdout is not None:
            stdout.close()
