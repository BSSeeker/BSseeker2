import gzip
import os

import datetime
import re
import shlex
import shutil
import string
import subprocess
import types
from itertools import izip

import marshal
import sys



_rc_trans = string.maketrans('ACGT', 'TGCA')
def reverse_compl_seq(strseq):
    return strseq.translate(_rc_trans)[::-1]



def show_version() :
    print ""
    print "     BS-Seeker2 v2.0.2 - April 20, 2013     "
    print ""

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

supported_aligners = [
                      BOWTIE,
                      BOWTIE2,
                      SOAP
                    ]

aligner_options_prefixes = { BOWTIE : '--bt-',
                             BOWTIE2 : '--bt2-',
                             SOAP   : '--soap-' }
aligner_path = dict((aligner, os.path.expanduser(find_location(aligner) or default_path)) for aligner, default_path in
                                                                                                   [(BOWTIE,'~/bowtie-0.12.7/'),
                                                                                                    (BOWTIE2, '~/bowtie-0.12.7/'),
                                                                                                    (SOAP, '~/soap2.21release/')
                                                                                                    ])


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
    input = (gzip.open if filename.endswith('.gz') else open)(filename, 'r')
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
    output.close()
    yield output_fname

    input.close()

def read_fasta(fasta_file):
    """ Iterates over all sequences in a fasta file. One at a time, without reading the whole file into the main memory.
    """

    input = (gzip.open if fasta_file.endswith('.gz') else open)(fasta_file)
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
        print "\n[Error]:\n\t Cannot fine file: %s.data" % filename
        exit(-1)

    obj = marshal.load(input)
    input.close()
    return obj   



def run_in_parallel(commands):

    commands = [(cmd[0], open(cmd[1], 'w')) if type(cmd) is tuple else (cmd, None) for cmd in commands]

    logm('Starting commands:\n' + '\n'.join([cmd for cmd, stdout in commands]))
    for i, proc in enumerate([subprocess.Popen(args = shlex.split(cmd), stdout = stdout) for cmd, stdout in commands]):
        return_code = proc.wait()
        logm('Finished: ' + commands[i][0])
        if return_code != 0:
            error('%s \nexited with an error code: %d. Please, check the log files.' % (commands[i], return_code))
    for _, stdout in commands:
        if stdout is not None:
            stdout.close()
