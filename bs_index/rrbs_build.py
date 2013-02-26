import os
from bs_utils.utils import *


FWD_MAPPABLE_REGIONS = lambda chrom_id: chrom_id+'.fwd_mappable_regions'
REV_MAPPABLE_REGIONS = lambda chrom_id: chrom_id+'.rev_mappable_regions'



def rrbs_build(fasta_file, build_command, ref_path, low_bound, up_bound, aligner):
    # ref_path is a string that containts the directory where the reference genomes are stored with
    # the input fasta filename appended
    ref_path = os.path.join(ref_path,
        os.path.split(fasta_file)[1] + '_rrbs_%d_%d' % (low_bound, up_bound) +'_' + aligner)

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

    for chrom_id, chrom_seq in read_fasta(fasta_file):
        total_chromosomes += 1
        refd[chrom_id] = len(chrom_seq)

        fwd_chr_regions = {}
        rev_chr_regions = {}

        L = len(chrom_seq)
        CCGG_sites = []
        CCGG_CCGG = []

        #-- collect all "CCGG sites ---------------------------------
        i = 1
        while i <= L - 4:
            if chrom_seq[i : i + 4] == "CCGG":
                CCGG_sites.append(i)
            i += 1

        #-- find "CCGG" pairs that are within the length of fragment ---------------------------------
        for j in xrange(len(CCGG_sites) - 1):
            DD = (CCGG_sites[j+1] - CCGG_sites[j]) - 4 # NOT including both CCGG
            if  low_bound <= DD <= up_bound:
                CCGG_CCGG.append([CCGG_sites[j], CCGG_sites[j+1] + 3]) # leftmost <--> rightmost
                mapable_seq = chrom_seq[CCGG_sites[j] : CCGG_sites[j+1] + 4]
                no_mappable_region += 1

                fwd_chr_regions[str(CCGG_sites[j])] = [CCGG_sites[j+1] + 3, no_mappable_region]
                rev_chr_regions[str(CCGG_sites[j+1] + 3)] = [CCGG_sites[j], no_mappable_region]

                # start_position, end_position, serial, sequence
                mappable_regions_output_file.write("%s\t%d\t%d\t%d\t%s\n"%(chrom_id, no_mappable_region, CCGG_sites[j], CCGG_sites[j+1]+3, mapable_seq))

        serialize(fwd_chr_regions, ref_p(FWD_MAPPABLE_REGIONS(chrom_id)))
        serialize(rev_chr_regions, ref_p(REV_MAPPABLE_REGIONS(chrom_id)))

        #-----------------------------------
        _map_seq = []
        mappable_length = 0
        unmappable_length = 0
        m = 0
        mark = False
        while m < L:
            if len(CCGG_CCGG) > 0:
                pair = CCGG_CCGG[0]
                p1 = pair[0]
                p2 = pair[1]
                if p1 <= m < p2 + 1:
                    _map_seq.append(chrom_seq[m])
                    mappable_length+=1
                    mark = True
                else:
                    if not mark:
                        _map_seq.append("-")
                        unmappable_length+=1
                    else: # the last eligible base
                        _ = CCGG_CCGG.pop(0)
                        if len(CCGG_CCGG)>0:
                            pair = CCGG_CCGG[0]
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

#
#
#
#    ###################################
#
#    genome, refd = fasta2dict(fasta_file)
#    serialize(refd, os.path.join(ref_path,"refname.data"))
#
#    all_base = sum(len(genome[key]) for key in genome)
#
#
#    logm("--- In total %d reference seqs ==> %d bp" % (len(genome),all_base) )
#
#
#    #--- positive strand -------------------------------------------------------------
#    FW_genome = {}
#    mappable_regions   = {}
#    mappable_regions_output_file = open(os.path.join(ref_path, "RRBS_mapable_regions.txt"),"w")
#
#    all_L = 0
#    all_mappable_length = 0
#    all_unmappable_length = 0
#
#    no_mappable_region = 0
#    for chr in sorted(genome.keys()):
#        chr_regions = []
#        seq = genome[chr]
#        L = len(seq)
#        CCGG_sites = []
#        CCGG_CCGG = []
#
#        #-- collect all "CCGG sites ---------------------------------
#        i = 1
#        while i <= L - 4:
#            if seq[i : i + 4] == "CCGG":
#                CCGG_sites.append(i)
#            i += 1
#
#        #-- find "CCGG" pairs that are within the length of fragment ---------------------------------
#        for j in xrange(len(CCGG_sites)-1):
#            DD = (CCGG_sites[j+1]-CCGG_sites[j])-4 # NOT including both CCGG
#            if  low_bound <= DD <= up_bound:
#                CCGG_CCGG.append([CCGG_sites[j], CCGG_sites[j+1]+3]) # leftmost <--> rightmost
#                mapable_seq = seq[CCGG_sites[j]:CCGG_sites[j+1]+4]
#                no_mappable_region += 1
#
#
#                chr_regions.append([CCGG_sites[j], CCGG_sites[j+1]+3, no_mappable_region, mapable_seq])
#                # start_position, end_position, serial, sequence
#                mappable_regions_output_file.write("%s\t%d\t%d\t%d\t%s\n"%(chr, no_mappable_region, CCGG_sites[j], CCGG_sites[j+1]+3, mapable_seq))
#
#        mappable_regions[chr] = chr_regions
#
#        #-----------------------------------
#        _map_seq = []
#        mappable_length = 0
#        unmappable_length = 0
#        m = 0
#        mark = False
#        while m < L:
#            if len(CCGG_CCGG) > 0:
#                pair = CCGG_CCGG[0]
#                p1 = pair[0]
#                p2 = pair[1]
#                if p1 <= m < p2 + 1:
#                    _map_seq.append(seq[m])
#                    mappable_length+=1
#                    mark = True
#                else:
#                    if not mark:
#                        _map_seq.append("-")
#                        unmappable_length+=1
#                    else: # the last eligible base
#                        _ = CCGG_CCGG.pop(0)
#                        if len(CCGG_CCGG)>0:
#                            pair = CCGG_CCGG[0]
#                            p1 = pair[0]
#                            p2 = pair[1]
#                            if  p1 <= m < p2 + 1:
#                                _map_seq.append(seq[m])
#                                mappable_length += 1
#                                mark = True
#                            else:
#                                _map_seq.append("-")
#                                unmappable_length += 1
#                                mark = False
#                        else:
#                            _map_seq.append("-")
#                            unmappable_length+=1
#                            mark = False
#            else:
#                if not mark:
#                    _map_seq.append("-")
#                    unmappable_length+=1
#                else:
#                    _map_seq.append("-")
#                    unmappable_length+=1
#                    mark = False
#
#            m+=1
#
#        del genome[chr]
#        #-----------------------------------
#
#        FW_genome[chr] = ''.join(_map_seq)
#
#        #-----------------------------------
#        logm("# %s : all %d : %d (unmapable -) %d (mapable) (%1.5f)"%(chr,
#                                                                               L,
#                                                                               unmappable_length,
#                                                                               mappable_length,
#                                                                               float(mappable_length)/L) )
#        all_L += L
#        all_mappable_length += mappable_length
#        all_unmappable_length += unmappable_length
#
#        elapsed(chr)
#
#    logm("# total %d chromosome seqs ==> %d : %d (unmapable -) %d (mapable) (%1.5f)" %(len(genome.keys()),
#                                                                                                all_L,
#                                                                                                all_unmappable_length,
#                                                                                                all_mappable_length,
#                                                                                                float(all_mappable_length)/all_L) )
#    logm("#       %d eligible fragments" % no_mappable_region )
#
#    elapsed('Cutting mappable regions')
#    serialize(FW_genome, os.path.join(ref_path, "ref.data"))
#    serialize(mappable_regions, os.path.join(ref_path, "RRBS_mapable_regions.data"))
#
##    outf.close()
#    mappable_regions_output_file.close()
#    elapsed('Store mappable regions and genome')
#
#
#    # Part 2
#    #----------------------------------------------------------------
#    logm("\n")
#    logm("----------------         Pre-processing mapable genome         (2-2) ----------------" )
#
#    FW_lst = sorted(FW_genome.iterkeys())
#
#    elapsed('Storing forward genome')
#
#
#    #---------------- 4 converted fasta -------------------------------------------
#
#    outf=open(os.path.join(ref_path,'W_C2T.fa'),'w')
#    for header in FW_lst:
#        outf.write('>%s\n%s\n' % (header, FW_genome[header].replace("C", "T")))
#    outf.close()
#    elapsed('end 4-1')
#
#    outf=open(os.path.join(ref_path,'W_G2A.fa'),'w')
#    for header in FW_lst:
#        outf.write('>%s\n%s\n' % (header, FW_genome[header].replace("G", "A")))
#    outf.close()
#    elapsed('end 4-3')
#
#
#    #---------------- Reverse complement (Crick strand) ----------------------------
#
#    # reverse complement in place to save memory
#    RC_genome = FW_genome
#    for key in RC_genome:
#        RC_genome[key] = reverse_compl_seq(RC_genome[key])
#    RC_lst=sorted(RC_genome.iterkeys())
#
#    outf=open(os.path.join(ref_path,'C_C2T.fa'),'w')
#    for header in RC_lst:
#        outf.write('>%s\n%s\n' % (header, RC_genome[header].replace("C", "T")))
#    outf.close()
#    elapsed('end 4-2')
#
#
#    outf=open(os.path.join(ref_path,'C_G2A.fa'),'w')
#    for header in RC_lst:
#        outf.write('>%s\n%s\n' % (header, RC_genome[header].replace("G", "A")))
#    outf.close()
#    elapsed('end 4-4')
#
#    #---------------- bowtie-build -------------------------------------------
#
#    # append ref_path to all elements of to_bowtie
#    to_bowtie = map(lambda f: os.path.join(ref_path, f), ['W_C2T', 'W_G2A', 'C_C2T', 'C_G2A'])
#
#    run_in_parallel([(build_command % { 'fname' : fname }, fname + '.log') for fname in to_bowtie])
#
#    elapsed('Index building')
#    # delete all fasta files
#    delete_files(f+'.fa' for f in to_bowtie)
#
##    gzip.wait()
#    elapsed('END')
