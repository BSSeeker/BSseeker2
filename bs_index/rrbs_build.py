import os

from bs_utils.utils import *


FWD_MAPPABLE_REGIONS = lambda chrom_id: chrom_id+'.fwd_mappable_regions'
REV_MAPPABLE_REGIONS = lambda chrom_id: chrom_id+'.rev_mappable_regions'

# example: T-CGA

def rrbs_build(fasta_file, build_command, ref_path, low_bound, up_bound, aligner, cut_format="C-CGG"):
    # ref_path is a string that contains the directory where the reference genomes are stored with
    # the input fasta filename appended

    cut_format = cut_format.upper()
    if cut_format == "C-CGG" :
        ref_path = os.path.join(ref_path,
            os.path.split(fasta_file)[1] + '_rrbs_%d_%d' % (low_bound, up_bound) +'_' + aligner)
    else :
        ref_path = os.path.join(ref_path,
            os.path.split(fasta_file)[1] + '_rrbs_%s_%d_%d' % (cut_format, low_bound, up_bound) +'_' + aligner)

    ref_p = lambda filename: os.path.join(ref_path, filename)


    clear_dir(ref_path)
    open_log(ref_p( 'log'))

    refd = {}
    w_c2t = open(ref_p('W_C2T.fa'),'w')
    c_c2t = open(ref_p('C_C2T.fa'),'w')

    w_g2a = open(ref_p('W_G2A.fa'),'w')
    c_g2a = open(ref_p('C_G2A.fa'),'w')

    mappable_regions_output_file = open(ref_p("RRBS_mapable_regions.txt"),"w")

    all_L = 0
    all_mappable_length = 0
    all_unmappable_length = 0

    no_mappable_region = 0

    total_chromosomes = 0

    cut_context = re.sub("-", "", cut_format)
    cut5_context = re.match( r'(.*)\-(.*)', cut_format).group(1) # 5' end
    cut3_context = re.match( r'(.*)\-(.*)', cut_format).group(2) # 3' end
    cut_len = len(cut_context)
    cut5_len = len(cut5_context)
    cut3_len = len(cut3_context)

    for chrom_id, chrom_seq in read_fasta(fasta_file):
        total_chromosomes += 1
        refd[chrom_id] = len(chrom_seq)

        fwd_chr_regions = {}
        rev_chr_regions = {}

        L = len(chrom_seq)
        XXXX_sites = []
        XXXX_XXXX = []

        #-- collect all "CCGG sites ---------------------------------
        i = 1
        while i <= L - cut_len:
            if chrom_seq[i : i + cut_len] == cut_context:
                XXXX_sites.append(i)
            i += 1

        #-- find "CCGG" pairs that are within the length of fragment ----
        for j in xrange(len(XXXX_sites) - 1):
            DD = (XXXX_sites[j+1] - XXXX_sites[j]) - cut_len # NOT including both CCGG; DD: fragment length
            if  low_bound <= DD <= up_bound:
#                XXXX_XXXX.append([XXXX_sites[j], XXXX_sites[j+1] + 3]) # leftmost <--> rightmost
                XXXX_XXXX.append([XXXX_sites[j], XXXX_sites[j+1] + cut_len - 1]) # leftmost <--> rightmost
                mapable_seq = chrom_seq[XXXX_sites[j] : XXXX_sites[j+1] + cut_len]
                no_mappable_region += 1

                fwd_chr_regions[str(XXXX_sites[j])] = [XXXX_sites[j+1] + cut_len - 1, no_mappable_region]
                rev_chr_regions[str(XXXX_sites[j+1] + cut_len -1)] = [XXXX_sites[j], no_mappable_region]

                # start_position, end_position, serial, sequence
                mappable_regions_output_file.write("%s\t%d\t%d\t%d\t%s\n"%(chrom_id, no_mappable_region,
                                            XXXX_sites[j], XXXX_sites[j+1] + cut_len - 1, mapable_seq))
        # storing region information to file
        # format: A[left_CCGG_pos]=[right_CCGG_pos, number_of_mappable_region]
        serialize(fwd_chr_regions, ref_p(FWD_MAPPABLE_REGIONS(chrom_id)))
        serialize(rev_chr_regions, ref_p(REV_MAPPABLE_REGIONS(chrom_id)))

        #-----------------------------------
        # mask the genome
        _map_seq = []
        mappable_length = 0
        unmappable_length = 0
        m = 0
        mark = False
        while m < L:
            if len(XXXX_XXXX) > 0:
                pair = XXXX_XXXX[0]
                p1 = pair[0]  # left end of fragment
                p2 = pair[1]  # right end of fragment
                if p1 <= m < p2 + 1:
                    _map_seq.append(chrom_seq[m])
                    mappable_length+=1
                    mark = True
                else:
                    if not mark:
                        _map_seq.append("-")
                        unmappable_length+=1
                    else: # the last eligible base
                        _ = XXXX_XXXX.pop(0)
                        if len(XXXX_XXXX)>0:
                            pair = XXXX_XXXX[0]
                            p1 = pair[0]
                            p2 = pair[1]
                            if  p1 <= m < p2 + 1:
                                _map_seq.append(chrom_seq[m])
                                mappable_length += 1
                                mark = True
                            else:
                                _map_seq.append("-")
                                unmappable_length += 1
                                mark = False
                        else:
                            _map_seq.append("-")
                            unmappable_length+=1
                            mark = False
            else:
                if not mark:
                    _map_seq.append("-")
                    unmappable_length+=1
                else:
                    _map_seq.append("-")
                    unmappable_length+=1
                    mark = False

            m+=1

        #-----------------------------------

        chrom_seq = ''.join(_map_seq)
        serialize(chrom_seq, ref_p(chrom_id))

        w_c2t.write('>%s\n%s\n' % (chrom_id, chrom_seq.replace("C","T")))
        w_g2a.write('>%s\n%s\n' % (chrom_id, chrom_seq.replace("G","A")))

        chrom_seq = reverse_compl_seq(chrom_seq)

        c_c2t.write('>%s\n%s\n' % (chrom_id, chrom_seq.replace("C","T")))
        c_g2a.write('>%s\n%s\n' % (chrom_id, chrom_seq.replace("G","A")))

        #-----------------------------------
        logm("# %s : all %d : %d (unmapable -) %d (mapable) (%1.5f)"%(chrom_id,
                                                                      L,
                                                                      unmappable_length,
                                                                      mappable_length,
                                                                      float(mappable_length)/L) )
        all_L += L
        all_mappable_length += mappable_length
        all_unmappable_length += unmappable_length

        elapsed('Finished initial preprocessing of ' + chrom_id)


    for outf in [w_c2t, c_c2t, w_g2a, c_g2a]:
        outf.close()


    logm("# total %d chromosome seqs ==> %d : %d (unmapable -) %d (mapable) (%1.5f)" %(total_chromosomes,
                                                                                       all_L,
                                                                                       all_unmappable_length,
                                                                                       all_mappable_length,
                                                                                       float(all_mappable_length)/all_L) )
    logm("#       %d eligible fragments" % no_mappable_region )

    serialize(refd, ref_p("refname"))

    mappable_regions_output_file.close()
    elapsed('Storing mappable regions and genome')

    #---------------- bowtie-build -------------------------------------------

    # append ref_path to all elements of to_bowtie
    to_bowtie = map(lambda f: os.path.join(ref_path, f), ['W_C2T', 'W_G2A', 'C_C2T', 'C_G2A'])

    run_in_parallel([(build_command % { 'fname' : fname }, fname + '.log') for fname in to_bowtie])

    elapsed('Index building')
    # delete all fasta files
    delete_files(f+'.fa' for f in to_bowtie)

    elapsed('END')

