import json
import sys

__author__ = 'pf'
import random
import string

_rc_trans = string.maketrans('ACGT', 'TGCA')
def reverse_compl_seq(strseq):
    return strseq.translate(_rc_trans)[::-1]


def generate_reads(out, cell_type, seq, nreads):
    nucs = 'ACTG'
    read_nucs = ['']*READ_LENGTH
    for i in xrange(nreads):
        
        read_seq = None
        while read_seq is None or 'N' in read_seq:
            start = random.randint(0, len(seq) - READ_LENGTH)
            end = start + READ_LENGTH
            read_seq = seq[start:end]
        
        strand = '+'
#        if random.random() < 0.5:
#            read_seq = reverse_compl_seq(read_seq)
#            strand = '-'
#
        mutations = 0
        for i, n in enumerate(read_seq):
            if random.random() < INDEL:
                if random.random() < 0.5:
                    read_nucs[i] = ''
                else:
                    read_nucs[i] = random.choice(nucs) + n

            elif random.random() > ERROR_RATE:
                read_nucs[i] = n
            else:
                read_nucs[i] = random.choice([k for k in nucs if k != n])
                mutations += 1


        read_seq = ''.join(read_nucs)
        if random.random() < ADAPTER_CONTAMINATION:
            read_seq = 'TCTGT' + read_seq
#        read_seq = ''.join(n if random.random() > ERROR_RATE else (random.choice([k for k in nucs if k != n])) for n in read_seq)
            
        
        out.write('>read_%s_%d_%s_%d_%d\n%s\n' % (cell_type, i, strand, start, end, read_seq))
        
        


if __name__ == '__main__':

    METHYLATION_RATE = 0.3
    COVERAGE = 6
    CELL_TYPES = 2
    READ_LENGTH = 75
    ERROR_RATE = 0.05

    INDEL = 0.01

    ADAPTER_CONTAMINATION = 0.1

#    inp_file = 'chr22.fa'
    inp_file = sys.argv[1]
#    inp_file = 'test_chr22.fa'

    genome = ''.join([l.strip() for l in list(open(inp_file)) if l[0] != '>']).upper()

    nreads = COVERAGE*len(genome)/READ_LENGTH

#    out = open('chr22_reads.fa','w')

    out = open(inp_file+'.reads.fa','w')

    print "reading genome done"
    methylation = {}

    methylated_nucs = ['']*len(genome)
    def methylate(seq, is_fwd_strand):

        for i in xrange(len(seq)):
            if seq[i] != 'C':
                methylated_nucs[i]= seq[i]
            else:
                pos = i if is_fwd_strand else len(seq) - i - 1

                if pos not in methylation:
                    methylation[pos] = {'C' : 0, 'T' : 0}

                if random.random() < METHYLATION_RATE:
                    methylated_nucs[i] = 'C'
                    methylation[pos]['C'] += 1
                else:
                    methylated_nucs[i] = 'T'
                    methylation[pos]['T'] += 1
        print 'methylating strand done'
        return ''.join(methylated_nucs)

    genome_rc = reverse_compl_seq(genome)

    for cell_type in xrange(CELL_TYPES):
        print cell_type

        generate_reads( out,
                        'fwd_%d' % cell_type,
                        methylate(genome, True),
                        nreads/2 )
                        
        generate_reads( out,
                        'rev_%d' % cell_type,
                        methylate(genome_rc, False),
                        nreads/2 )

    for i in methylation:
        methylation[i]['level'] = float(methylation[i]['C'])/(methylation[i]['C'] + methylation[i]['T']) if (methylation[i]['C'] + methylation[i]['T']) > 0 else 0
    json.dump(methylation.items(), open(inp_file + '_methylation.json','w'))
        
    
