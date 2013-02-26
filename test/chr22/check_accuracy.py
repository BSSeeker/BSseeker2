import sys
import json
import math

import pysam

def bs_seeker(outf):
    for a in pysam.Samfile(outf, 'rb'):
        _, original_strand, ct, rno, strand, start, end  = a.qname.split('_')
        yield a.qname, original_strand, strand, start, end, a.pos

def bismark(outf):

    for l in open(outf):
        if l[0] == '@':
            continue
        buf = l.split()
        _, original_strand, ct, rno, strand, start, end  = buf[0].split('_')
        yield  buf[0], original_strand, strand, start, end, buf[3]



def check_cgmap(cgmap_file):
    from scipy.stats.stats import pearsonr
#    meth_levels = json.load(open('chr22_methylation.json'))
    meth_levels = json.load(open('chr22.fa_methylation.json'))

    print 'reading meth levels done'

    parse_ml = lambda a: (int(a[2]) - 1, float(a[5]))
    cgmap = dict(parse_ml(l.split()) for l in cgmap_file)

    print 'reading cgmap done'

    error = 0
    e2 = 0
    h = dict(meth_levels)
    for pos in cgmap:
        if pos not in h:
            print pos

    for pos, info in meth_levels:

        if pos not in cgmap:
            cgmap[pos] = 0

        error += (info['level'] - cgmap[pos])**2
        e2 += abs(info['level'] - cgmap[pos])

    print len(cgmap), len(meth_levels)
    print 'RMSE:', math.sqrt(error/len(meth_levels)), e2/len(meth_levels)

    print 'pearson r:', pearsonr([v for k,v in sorted(cgmap.iteritems())],
                                 [v for k,v in sorted((pos, info['level']) for pos, info in meth_levels)])



def rmapbs(outf):
    pass


if __name__ == '__main__':


    if len(sys.argv) != 3:
        print "usage: %s bs_seeker|bismark|rmapbs|cgmap output-file" % __file__
        exit(1)

    total = 0
    tp = 0
    if sys.argv[1] == 'bs_seeker':
        parse_output = bs_seeker
    elif sys.argv[1] == 'bismark':
        parse_output = bismark
    elif sys.argv[1] == 'rmapbs':
        parse_output = rmapbs
    elif sys.argv[1] == 'cgmap':
        check_cgmap(open(sys.argv[2]))
        exit(0)
    else:
        print "usage: %s bs_seeker|bismark|rmapbs|cgmap output-file" % __file__
        exit(1)

    chr22_len = 51304566

    for rid, original_strand, strand, start, end, mapped_pos in parse_output(sys.argv[2]):
        start, end, mapped_pos = map(int, [start, end, mapped_pos])
        if original_strand == 'rev':
            start = chr22_len - start
            end = chr22_len - end
            start, end = end, start
        if abs(start - mapped_pos) < 20:
            tp += 1
#            print rid, start, end, mapped_pos
#        else:
#            print rid, start, end, mapped_pos
        total += 1
    print 'total:', total, 'tp:', tp, '(%.2lf)' % (100*float(tp)/total)

