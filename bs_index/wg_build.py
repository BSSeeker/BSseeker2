from bs_utils.utils import *


def wg_build(fasta_file, build_command, ref_path, aligner):

    # ref_path is a string that contains the directory where the reference genomes are stored with
    # the input Fasta filename appended
    ref_path = os.path.join(ref_path,
                            os.path.split(fasta_file)[1] + '_'+aligner)

    clear_dir(ref_path)
    #---------------------------------------------------------------
    # 1. First get the complementary genome (also do the reverse)
    # 2. Then do CT and GA conversions
    #---------------------------------------------------------------

    open_log(os.path.join(ref_path, 'log'))
    refd = {}
    w_c2t = open(os.path.join(ref_path, 'W_C2T.fa'),'w')
    c_c2t = open(os.path.join(ref_path, 'C_C2T.fa'),'w')

    w_g2a = open(os.path.join(ref_path, 'W_G2A.fa'),'w')
    c_g2a = open(os.path.join(ref_path, 'C_G2A.fa'),'w')
    for chrom_id, chrom_seq in read_fasta(fasta_file):
        serialize(chrom_seq, os.path.join(ref_path, chrom_id))
        refd[chrom_id] = len(chrom_seq)

        w_c2t.write('>%s\n%s\n' % (chrom_id, chrom_seq.replace("C","T")))
        w_g2a.write('>%s\n%s\n' % (chrom_id, chrom_seq.replace("G","A")))

        chrom_seq = reverse_compl_seq(chrom_seq)

        c_c2t.write('>%s\n%s\n' % (chrom_id, chrom_seq.replace("C","T")))
        c_g2a.write('>%s\n%s\n' % (chrom_id, chrom_seq.replace("G","A")))

        elapsed('Preprocessing '+chrom_id)

    for outf in [w_c2t, c_c2t, w_g2a, c_g2a]:
        outf.close()

    serialize(refd, os.path.join(ref_path,"refname"))
    elapsed('Genome preprocessing')
    # append ref_path to all elements of to_bowtie
    to_bowtie = map(lambda f: os.path.join(ref_path, f), ['W_C2T', 'W_G2A', 'C_C2T', 'C_G2A'])

    # start bowtie-build for all converted genomes and wait for the processes to finish

    #run_in_parallel([(build_command % { 'fname' : fname }, fname+'.log') for fname in to_bowtie])
    # --------------- revised on 2018-01-30
    for fname in to_bowtie :
        run_in_parallel([(build_command % { 'fname' : fname }, fname+'.log')])
    # ---------------
    # delete fasta files of converted genomes
    if aligner != "rmap" :
        delete_files(f+'.fa' for f in to_bowtie)

    elapsed('Done')

