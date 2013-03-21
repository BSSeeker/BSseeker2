__author__ = 'pf'

from bs_utils.utils import *
import re


BAM_MATCH = 0
BAM_INS = 1
BAM_DEL = 2
BAM_SOFTCLIP = 4

CIGAR_OPS = {'M' : BAM_MATCH, 'I' : BAM_INS, 'D' : BAM_DEL, 'S' : BAM_SOFTCLIP}


def N_MIS(r,g):
    mismatches = 0
    if len(r)==len(g):
        for i in xrange(len(r)):
            if r[i] != g[i] and r[i] != "N" and g[i] != "N" and not(r[i] == 'T' and g[i] == 'C'):
                mismatches += 1
    return mismatches


#----------------------------------------------------------------

def next_nuc(seq, pos, n):
    """ Returns the nucleotide that is n places from pos in seq. Skips gap symbols.
    """
    i = pos + 1
    while i < len(seq):
        if seq[i] != '-':
            n -= 1
            if n == 0: break
        i += 1
    return seq[i]



def methy_seq(read, genome):
    H = ['A', 'C', 'T']
    m_seq = []
    xx = "-"
    for i in xrange(len(read)):

        if genome[i] == '-':
            continue

        elif read[i] != 'C' and read[i] != 'T':
            xx = "-"

        elif read[i] == "T" and genome[i] == "C": #(unmethylated):
            nn1 = next_nuc(genome, i, 1)
            if nn1 == "G":
                xx = "x"
            elif nn1 in H :
                nn2 = next_nuc(genome, i, 2)
                if nn2 == "G":
                    xx = "y"
                elif nn2 in H :
                    xx = "z"

        elif read[i] == "C" and genome[i] == "C": #(methylated):
            nn1 = next_nuc(genome, i, 1)

            if nn1 == "G":
                xx = "X"
            elif nn1 in H :
                nn2 = next_nuc(genome, i, 2)

                if nn2 == "G":
                    xx = "Y"
                elif nn2 in H:
                    xx = "Z"
        else:
            xx = "-"
        m_seq.append(xx)

    return ''.join(m_seq)

def mcounts(mseq, mlst, ulst):
    out_mlst=[mlst[0]+mseq.count("X"), mlst[1]+mseq.count("Y"), mlst[2]+mseq.count("Z")]
    out_ulst=[ulst[0]+mseq.count("x"), ulst[1]+mseq.count("y"), ulst[2]+mseq.count("z")]
    return out_mlst, out_ulst



def process_aligner_output(filename, pair_end = False):

    #m = re.search(r'-('+'|'.join(supported_aligners) +')-TMP', filename)
    m = re.search(r'-('+'|'.join(supported_aligners) +')-.*TMP', filename)
    if m is None:
        error('The temporary folder path should contain the name of one of the supported aligners: ' + filename)

    format = m.group(1)

    input = open(filename)

    QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL = range(11)
    def parse_SAM(line):
        buf = line.split()

        flag = int(buf[FLAG])

        # skip reads that are not mapped
        # skip reads that have probability of being non-unique higher than 1/10
        if flag & 0x4 : # or int(buf[MAPQ]) < 10:
            return None, None, None, None, None, None
        # print "format = ", format
        if format == BOWTIE:
            mismatches = int([buf[i][5:] for i in xrange(11, len(buf)) if buf[i][:5] == 'NM:i:'][0]) # get the edit distance
        # --- bug fixed ------
        elif format == BOWTIE2:
            if re.search(r'(.)*-e2e-TMP(.*)', filename) is None : # local model
                mismatches = 1-int([buf[i][5:] for i in xrange(11, len(buf)) if buf[i][:5] == 'AS:i:'][0])
                # print "====local=====\n"
                ## bowtie2 use AS tag (score) to evaluate the mapping. The higher, the better.
            else : # end-to-end model
                # print "end-to-end\n"
                mismatches = int([buf[i][5:] for i in xrange(11, len(buf)) if buf[i][:5] == 'XM:i:'][0])
        # --- Weilong ---------
        else:
            mismatches = 1-buf[MAPQ]
            # mismatches = 1/float(buf[MAPQ])
            ## downstream might round (0,1) to 0, so use integer instead
            ## fixed by Weilong

#        # add the soft clipped nucleotides to the number of mismatches
#        cigar_string = buf[CIGAR]
#        cigar = parse_cigar(cigar_string)
#        if cigar[0][0] == BAM_SOFTCLIP:
#            mismatches += cigar[0][1]
#        if cigar[-1][0] == BAM_SOFTCLIP:
#            mismatches += cigar[-1][1]

        return (buf[QNAME], # read ID
                buf[RNAME], # reference ID
                int(buf[POS]) - 1, # position, 0 based (SAM is 1 based)
                mismatches,    # number of mismatches
                parse_cigar(buf[CIGAR]), # the parsed cigar string
                flag & 0x40 # true if it is the first mate in a pair, false if it is the second mate
                )

    SOAP_QNAME, SOAP_SEQ, SOAP_QUAL, SOAP_NHITS, SOAP_AB, SOAP_LEN, SOAP_STRAND, SOAP_CHR, SOAP_LOCATION, SOAP_MISMATCHES = range(10)
    def parse_SOAP(line):
        buf = line.split()
        return (buf[SOAP_QNAME],
                buf[SOAP_CHR],
                int(buf[SOAP_LOCATION]) - 1,
                int(buf[SOAP_MISMATCHES]),
                buf[SOAP_AB],
                buf[SOAP_STRAND],
                parse_cigar(buf[SOAP_LEN]+'M')
            )


    if format == BOWTIE or format == BOWTIE2:
        if pair_end:
            for line in input:
                header1, chr1, location1, no_mismatch1, cigar1,        _ = parse_SAM(line)
                header2,    _, location2, no_mismatch2, cigar2, mate_no2 = parse_SAM(input.next())


                if header1 and header2:
                    # flip the location info if the second mate comes first in the alignment file
                    if mate_no2:
                        location1, location2 = location2, location1
                        cigar1, cigar2 = cigar2, cigar1


                    yield header1, chr1, no_mismatch1 + no_mismatch2, location1, cigar1, location2, cigar2
        else:
            for line in input:
                header, chr, location, no_mismatch, cigar, _ = parse_SAM(line)
                if header is not None:
                    yield header, chr, location, no_mismatch, cigar
    elif format == SOAP:
        if pair_end:
            for line in input:
                header1, chr1, location1, no_mismatch1, mate1, strand1, cigar1 = parse_SOAP(line)
                header2, _   , location2, no_mismatch2,     _, strand2, cigar2 = parse_SOAP(input.next())

                if mate1 == 'b':
                    location1, location2 = location2, location1
                    strand1, strand2 = strand2, strand1
                    ciga1, cigar2 = cigar2, cigar1


                if header1 and header2 and strand1 == '+' and strand2 == '-':
                    yield header1, chr1, no_mismatch1 + no_mismatch2, location1, cigar1, location2, cigar2

        else:
            for line in input:
                header, chr, location, no_mismatch, _, strand, cigar = parse_SOAP(line)
                if header and strand == '+':
                    yield header, chr, location, no_mismatch, cigar

    input.close()


def parse_cigar(cigar_string):
    i = 0
    prev_i = 0
    cigar = []
    while i < len(cigar_string):
        if cigar_string[i] in CIGAR_OPS:
            cigar.append((CIGAR_OPS[cigar_string[i]], int(cigar_string[prev_i:i])))
            prev_i = i + 1
        i += 1
    return cigar

def get_read_start_end_and_genome_length(cigar):
    r_start = cigar[0][1] if cigar[0][0] == BAM_SOFTCLIP else 0
    r_end = r_start
    g_len = 0
    for edit_op, count in cigar:
        if edit_op == BAM_MATCH:
            r_end += count
            g_len += count
        elif edit_op == BAM_INS:
            r_end += count
        elif edit_op == BAM_DEL:
            g_len += count
    return r_start, r_end, g_len # return the start and end in the read and the length of the genomic sequence



def cigar_to_alignment(cigar, read_seq, genome_seq):
    """ Reconstruct the pairwise alignment based on the CIGAR string and the two sequences
    """

    # reconstruct the alignment
    r_pos = cigar[0][1] if cigar[0][0] == BAM_SOFTCLIP else 0
    g_pos = 0
    r_aln = ''
    g_aln = ''
    for edit_op, count in cigar:
        if edit_op == BAM_MATCH:
            r_aln += read_seq[r_pos : r_pos + count]
            g_aln += genome_seq[g_pos : g_pos + count]
            r_pos += count
            g_pos += count
        elif edit_op == BAM_DEL:
            r_aln += '-'*count
            g_aln += genome_seq[g_pos : g_pos + count]
            g_pos += count
        elif edit_op == BAM_INS:
            r_aln += read_seq[r_pos : r_pos + count]
            g_aln += '-'*count
            r_pos += count

    return r_aln, g_aln




def get_genomic_sequence(genome, start, end, strand = '+'):
    if start > 1:
        prev = genome[start-2:start]
    elif start == 1:
        prev = 'N'+genome[0]
    else:
        prev = 'NN'

    if end < len(genome) - 1:
        next = genome[end: end + 2]
    elif end == len(genome) - 1:
        next = genome[end] + 'N'
    else:
        next = 'NN'
    origin_genome = genome[start:end]

    if strand == '-':
        # reverse complement everything if strand is '-'
        revc = reverse_compl_seq('%s%s%s' % (prev, origin_genome, next))
        prev, origin_genome, next = revc[:2], revc[2:-2], revc[-2:]

    return origin_genome, next, '%s_%s_%s' % (prev, origin_genome, next)
