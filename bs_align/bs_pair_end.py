import fileinput, os, random, math
from bs_utils.utils import *
from bs_align_utils import *
#----------------------------------------------------------------
def extract_mapping(ali_file):
    unique_hits = {}
    non_unique_hits = {}
    header0 = ""
    family = []
    for header, chr, no_mismatch, location1, cigar1, location2, cigar2 in process_aligner_output(ali_file, pair_end = True):
        #------------------------
        if header != header0:

            # --- output ----
            if len(family) == 1:
                unique_hits[header0] = family[0]
            elif len(family) > 1:
                min_lst = min(family, key = lambda x: x[0])
                max_lst = max(family, key = lambda x: x[0])

                if min_lst[0] < max_lst[0]:
                    unique_hits[header0] = min_lst
                else:
                    non_unique_hits[header0] = min_lst[0]

            header0 = header
            family = []
        family.append((no_mismatch, chr, location1, cigar1, location2, cigar2))

    if len(family) == 1:
        unique_hits[header0] = family[0]
    elif len(family) > 1:
        min_lst = min(family, key = lambda x: x[0])
        max_lst = max(family, key = lambda x: x[0])

        if min_lst[0] < max_lst[0]:
            unique_hits[header0] = min_lst
        else:
            non_unique_hits[header0] = min_lst[0]

    return unique_hits, non_unique_hits

def _extract_mapping(ali_file):
    U = {}
    R = {}
    header0 = ""
    family = []
    for line in fileinput.input(ali_file):
        l = line.split()
        header = l[0][:-2]
        chr = str(l[1])
        location = int(l[2])
        #no_hits=int(l[4])
        #-------- mismatches -----------
        if len(l) == 4:
            no_mismatch = 0
        elif len(l) == 5:
            no_mismatch = l[4].count(":")
        else:
            print l
        #------------------------
        if header != header0:
            #--------------------
            if header0 != "":
                # --- output ----
                if len(family) == 1:
                    U[header0] = family[0]
                else:
                    if family[0][0] < family[1][0]:
                        U[header0] = family[0]
                    elif family[1][0] < family[0][0]:
                        U[header0] = family[1]
                    else:
                        R[header0] = family[0][0]
                family=[]
                # ---------------
            header0 = header
            family = [[no_mismatch, chr, location]]
            member = 1
        elif header == header0:
            if member == 1:
                family[-1][0] += no_mismatch
                family[-1].append(location)
                member = 2
            elif member == 2:
                family.append([no_mismatch, chr, location])
                member = 1
        #------------------------------

    fileinput.close()
    return U, R


#----------------------------------------------------------------

def bs_pair_end(main_read_file_1,
                main_read_file_2,
                asktag,
                adapter_file,
                cut1,
                cut2, # add cut3 and cut4?
                no_small_lines,
                max_mismatch_no,
                aligner_command,
                db_path,
                tmp_path,
                outfile, XS_pct, XS_count, show_multiple_hit=False):


    #----------------------------------------------------------------
    adapter=""
    adapterA=""
    adapterB=""
    if adapter_file !="":
        try :
            adapter_inf=open(adapter_file,"r")
            if asktag=="N": #<--- directional library
                adapter=adapter_inf.readline()
                adapter_inf.close()
                adapter=adapter.rstrip("\n")
            elif asktag=="Y":#<--- un-directional library
                adapterA=adapter_inf.readline()
                adapterB=adapter_inf.readline()
                adapter_inf.close()
                adapterA=adapterA.rstrip("\n")
                adapterB=adapterB.rstrip("\n")
        except IOError :
            print "[Error] Cannot find adapter file : %s" % adapter_file
            exit(-1)


    #----------------------------------------------------------------

    logm("End 1 filename: %s"% main_read_file_1  )
    logm("End 2 filename: %s"% main_read_file_2  )
    logm("The first base (for mapping): %d"% cut1  )
    logm("The last base (for mapping): %d"% cut2  )

    logm("-------------------------------- " )
    logm("Un-directional library: %s" % asktag  )
    logm("Path for short reads aligner: %s"% aligner_command + '\n')
    logm("Reference genome library path: %s"% db_path  )
    logm("Number of mismatches allowed: %s"% max_mismatch_no  )

    if adapter_file !="":
        if asktag=="Y":
            logm("Adapters to be removed from 3' of the reads:" )
            logm("-- A: %s" % adapterA  )
            logm("-- B: %s" % adapterB  )
        elif asktag=="N":
            logm("Adapter to be removed from 3' of the reads:" )
            logm("-- %s" % adapter  )

    logm("-------------------------------- " )

    #----------------------------------------------------------------

    # helper method to join fname with tmp_path
    tmp_d = lambda fname: os.path.join(tmp_path, fname)

    db_d = lambda fname: os.path.join(db_path, fname)


    #----------------------------------------------------------------
    # splitting the 2 big read files

    input_fname1 = os.path.split(main_read_file_1)[1]
    input_fname2 = os.path.split(main_read_file_2)[1]

    # TODO: run these in parallel with a subprocess
    split_file(main_read_file_1, tmp_d(input_fname1)+'-E1-', no_small_lines)
    split_file(main_read_file_2, tmp_d(input_fname2)+'-E2-', no_small_lines)

    dirList=os.listdir(tmp_path)
    my_files = zip(sorted(filter(lambda fname: fname.startswith("%s-E1-" % input_fname1), dirList)),
                   sorted(filter(lambda fname: fname.startswith("%s-E2-" % input_fname2), dirList)))



    #---- Stats ------------------------------------------------------------
    all_raw_reads=0
    all_trimmed=0
    all_mapped=0
    all_mapped_passed=0

    all_unmapped=0

    numbers_premapped_lst=[0,0,0,0]
    numbers_mapped_lst=[0,0,0,0]

    mC_lst=[0,0,0]
    uC_lst=[0,0,0]

    no_my_files=0

    #----------------------------------------------------------------
    print "== Start mapping =="

    for read_file_1, read_file_2 in my_files:

        no_my_files+=1

        random_id=".tmp-"+str(random.randint(1000000,9999999))

        original_bs_reads_1 = {}
        original_bs_reads_2 = {}
        original_bs_reads_lst= [original_bs_reads_1, original_bs_reads_2]


        if asktag == "Y" :

            #----------------------------------------------------------------
            outfile_1FCT = tmp_d('Trimmed_FCT_1.fa'+random_id)
            outfile_1RCT = tmp_d('Trimmed_RCT_1.fa'+random_id)
            outfile_2FCT = tmp_d('Trimmed_FCT_2.fa'+random_id)
            outfile_2RCT = tmp_d('Trimmed_RCT_2.fa'+random_id)

            try :
                read_inf = open(tmp_d(read_file_1),"r")
            except IOError :
                print "[Error] Cannot open file : %s", tmp_d(read_file_1)
                exit(-1)
            oneline = read_inf.readline()
            l = oneline.split()
            input_format = ""

            if oneline[0]=="@":	# Illumina GAII FastQ (Lister et al Nature 2009)
                input_format="fastq"
                n_fastq=0
            elif len(l)==1 and oneline[0]!=">": 	# pure sequences
                input_format="seq"
            elif len(l)==11:	# Illumina GAII qseq file
                input_format="qseq"
            elif oneline[0]==">":	# fasta
                input_format="fasta"
                n_fasta=0
            read_inf.close()

            print "Detected data format: %s" % input_format

            #----------------------------------------------------------------
            read_file_list   = [read_file_1, read_file_2]
            outfile_FCT_list = [outfile_1FCT, outfile_2FCT]
            outfile_RCT_list = [outfile_1RCT, outfile_2RCT]
            n_list = [0, 0]


            for f in range(2):
                read_file = read_file_list[f]
                outf_FCT = open(outfile_FCT_list[f], 'w')
                outf_RCT = open(outfile_RCT_list[f], 'w')
                original_bs_reads = original_bs_reads_lst[f]
                n = n_list[f]

                seq_id = ""
                seq = ""
                seq_ready = "N"
                for line in fileinput.input(tmp_d(read_file)):
                    l=line.split()
                    if input_format=="seq":
                        n+=1
                        seq_id=str(n)
                        seq_id=seq_id.zfill(12)
                        seq=l[0]
                        seq_ready="Y"
                    elif input_format=="fastq":
                        m_fastq=math.fmod(n_fastq,4)
                        n_fastq+=1
                        seq_ready="N"
                        if m_fastq==0:
                            n+=1
                            seq_id=str(n)
                            seq_id=seq_id.zfill(12)
                            seq=""
                        elif m_fastq==1:
                            seq=l[0]
                            seq_ready="Y"
                        else:
                            seq=""
                    elif input_format=="qseq":
                        n+=1
                        seq_id=str(n)
                        seq_id=seq_id.zfill(12)
                        seq=l[8]
                        seq_ready="Y"
                    elif input_format=="fasta":
                        m_fasta=math.fmod(n_fasta,2)
                        n_fasta+=1
                        seq_ready="N"
                        if m_fasta==0:
                            n+=1
                            seq_id=l[0][1:]
                            seq_id=seq_id.zfill(17)
                            seq=""
                        elif m_fasta==1:
                            seq=l[0]
                            seq_ready="Y"
                        else:
                            seq=""
                    #----------------------------------------------------------------
                    if seq_ready=="Y":
                        seq=seq[cut1-1:cut2] #------- selecting 0..52 from 1..72  -e 52
                        seq=seq.upper()
                        seq=seq.replace(".","N")

                        #--striping BS adapter from 3' read --------------------------------------------------------------
                        if (adapterA !="") and (adapterB !=""):
                            signature=adapterA[:6]
                            if signature in seq:
                                signature_pos=seq.index(signature)
                                if seq[signature_pos:] in adapterA:
                                    seq=seq[:signature_pos]#+"".join(["N" for x in range(len(seq)-len(signature_pos))])
                                    all_trimmed+=1
                            else:
                                signature=adapterB[:6]
                                if signature in seq:
                                    #print seq_id,seq,signature;
                                    signature_pos=seq.index(signature)
                                    if seq[signature_pos:] in adapterB:
                                        seq=seq[:signature_pos]#+"".join(["N" for x in range(len(seq)-len(signature_pos))])
                                        all_trimmed+=1

                        if len(seq) <= 4:
                            seq = "N" * (cut2-cut1+1)
                        #---------  trimmed_raw_BS_read  ------------------
                        original_bs_reads[seq_id] = seq

                        #---------  FW_C2T  ------------------
                        outf_FCT.write('>%s\n%s\n' % (seq_id, seq.replace("C","T")))
                        #---------  RC_G2A  ------------------
                        outf_RCT.write('>%s\n%s\n' % (seq_id, seq.replace("G","A")))

                n_list[f]=n

                outf_FCT.close()
                outf_RCT.close()

                fileinput.close()


            #print "All input end 1: %d , end 2: %d "%(n_list[0],n_list[1]);
            all_raw_reads+=n

            #--------------------------------------------------------------------------------
            # Bowtie mapping
            #--------------------------------------------------------------------------------
            WC2T_fr=tmp_d("W_C2T_fr_m"+max_mismatch_no+".mapping"+random_id)
            WC2T_rf=tmp_d("W_C2T_rf_m"+max_mismatch_no+".mapping"+random_id)
            CC2T_fr=tmp_d("C_C2T_fr_m"+max_mismatch_no+".mapping"+random_id)
            CC2T_rf=tmp_d("C_C2T_rf_m"+max_mismatch_no+".mapping"+random_id)

            run_in_parallel([aligner_command % {'reference_genome' : os.path.join(db_path,'W_C2T'),
                                                      'input_file_1' : outfile_1FCT,
                                                      'input_file_2' : outfile_2RCT,
                                                      'output_file' : WC2T_fr},

                             aligner_command % {'reference_genome' : os.path.join(db_path,'C_C2T'),
                                                      'input_file_1' : outfile_1FCT,
                                                      'input_file_2' : outfile_2RCT,
                                                      'output_file' : CC2T_fr},

                             aligner_command % {'reference_genome' : os.path.join(db_path,'W_C2T'),
                                                      'input_file_1' : outfile_2FCT,
                                                      'input_file_2' : outfile_1RCT,
                                                      'output_file' : WC2T_rf},

                             aligner_command % {'reference_genome' : os.path.join(db_path,'C_C2T'),
                                                      'input_file_1' : outfile_2FCT,
                                                      'input_file_2' : outfile_1RCT,
                                                      'output_file' : CC2T_rf}])


            delete_files(outfile_1FCT, outfile_2FCT, outfile_1RCT, outfile_2RCT)


            #--------------------------------------------------------------------------------
            # Post processing
            #--------------------------------------------------------------------------------


            FW_C2T_fr_U, FW_C2T_fr_R = extract_mapping(WC2T_fr)
            FW_C2T_rf_U, FW_C2T_rf_R = extract_mapping(WC2T_rf)
            RC_C2T_fr_U, RC_C2T_fr_R = extract_mapping(CC2T_fr)
            RC_C2T_rf_U, RC_C2T_rf_R = extract_mapping(CC2T_rf)

            delete_files(WC2T_fr, WC2T_rf, CC2T_fr, CC2T_rf)

            #----------------------------------------------------------------
            # get uniq-hit reads
            #----------------------------------------------------------------
            Union_set=set(FW_C2T_fr_U.iterkeys()) | set(FW_C2T_rf_U.iterkeys()) | set(RC_C2T_fr_U.iterkeys()) | set(RC_C2T_rf_U.iterkeys())

            Unique_FW_fr_C2T=set() # +
            Unique_FW_rf_C2T=set() # +
            Unique_RC_fr_C2T=set() # -
            Unique_RC_rf_C2T=set() # -
            Multiple_hits=set()


            for x in Union_set:
                _list=[]
                for d in [FW_C2T_fr_U, FW_C2T_rf_U, RC_C2T_fr_U, RC_C2T_rf_U]:
                    mis_lst=d.get(x,[99])
                    mis=int(mis_lst[0])
                    _list.append(mis)
                for d in [FW_C2T_fr_R, FW_C2T_rf_R, RC_C2T_fr_R, RC_C2T_rf_R]:
                    mis=d.get(x,99)
                    _list.append(mis)
                mini=min(_list)
                if _list.count(mini)==1:
                    mini_index=_list.index(mini)
                    if mini_index==0:
                        Unique_FW_fr_C2T.add(x)
                    elif mini_index==1:
                        Unique_FW_rf_C2T.add(x)
                    elif mini_index==2:
                        Unique_RC_fr_C2T.add(x)
                    elif mini_index==3:
                        Unique_RC_rf_C2T.add(x)
                    else :
                        Multiple_hits.add(x)
                else :
                    Multiple_hits.add(x)

            # write reads rejected by Multiple Hits to file
            if show_multiple_hit :
                outf_MH=open("Multiple_hit.fa",'w')
                for i in Multiple_hits :
                    outf_MH.write(">%s\n" % i)
                    outf_MH.write("%s\n" % original_bs_reads[i])
                outf_MH.close()

            del Union_set
            del FW_C2T_fr_R
            del FW_C2T_rf_R
            del RC_C2T_fr_R
            del RC_C2T_rf_R

            FW_C2T_fr_uniq_lst=[[FW_C2T_fr_U[u][1],u] for u in Unique_FW_fr_C2T]
            FW_C2T_rf_uniq_lst=[[FW_C2T_rf_U[u][1],u] for u in Unique_FW_rf_C2T]
            RC_C2T_fr_uniq_lst=[[RC_C2T_fr_U[u][1],u] for u in Unique_RC_fr_C2T]
            RC_C2T_rf_uniq_lst=[[RC_C2T_rf_U[u][1],u] for u in Unique_RC_rf_C2T]

            FW_C2T_fr_uniq_lst.sort()
            FW_C2T_rf_uniq_lst.sort()
            RC_C2T_fr_uniq_lst.sort()
            RC_C2T_rf_uniq_lst.sort()
            FW_C2T_fr_uniq_lst=[x[1] for x in FW_C2T_fr_uniq_lst]
            FW_C2T_rf_uniq_lst=[x[1] for x in FW_C2T_rf_uniq_lst]
            RC_C2T_fr_uniq_lst=[x[1] for x in RC_C2T_fr_uniq_lst]
            RC_C2T_rf_uniq_lst=[x[1] for x in RC_C2T_rf_uniq_lst]
            #----------------------------------------------------------------

            numbers_premapped_lst[0]+=len(Unique_FW_fr_C2T)
            numbers_premapped_lst[1]+=len(Unique_FW_rf_C2T)
            numbers_premapped_lst[2]+=len(Unique_RC_fr_C2T)
            numbers_premapped_lst[3]+=len(Unique_RC_rf_C2T)

            del Unique_FW_fr_C2T
            del Unique_FW_rf_C2T
            del Unique_RC_fr_C2T
            del Unique_RC_rf_C2T

            #----------------------------------------------------------------

            nn = 0
            for ali_unique_lst, ali_dic in [(FW_C2T_fr_uniq_lst,FW_C2T_fr_U),
                                            (FW_C2T_rf_uniq_lst,FW_C2T_rf_U),
                                            (RC_C2T_fr_uniq_lst,RC_C2T_fr_U),
                                            (RC_C2T_rf_uniq_lst,RC_C2T_rf_U)]:
                nn += 1

                mapped_chr0 = ""
                for header in ali_unique_lst:

                    _, mapped_chr, mapped_location_1, cigar1, mapped_location_2, cigar2 = ali_dic[header]

                    #-------------------------------------
                    if mapped_chr != mapped_chr0:
                        my_gseq=deserialize(db_d(mapped_chr))
                        chr_length=len(my_gseq)
                        mapped_chr0=mapped_chr
                    #-------------------------------------
                    if nn == 1 or nn == 3:
                        original_BS_1 = original_bs_reads_1[header]
                        original_BS_2 = reverse_compl_seq(original_bs_reads_2[header])
                    else:
                        original_BS_1 = original_bs_reads_2[header]
                        original_BS_2 = reverse_compl_seq(original_bs_reads_1[header])

                    r_start_1, r_end_1, g_len_1 = get_read_start_end_and_genome_length(cigar1)
                    r_start_2, r_end_2, g_len_2 = get_read_start_end_and_genome_length(cigar2)

                    all_mapped += 1

                    if nn == 1: 							# FW-RC mapped to + strand:

                        FR = "+FR"

#                        mapped_location_1 += 1
#                        origin_genome_long_1 = my_gseq[mapped_location_1 - 2 - 1 : mapped_location_1 + g_len_1 + 2 - 1]
#                        origin_genome_1 = origin_genome_long_1[2:-2]
                        mapped_strand_1 = "+"

#                        mapped_location_2 += 1
#                        origin_genome_long_2 = my_gseq[mapped_location_2 - 2 - 1 : mapped_location_2 + g_len_2 + 2 - 1]
#                        origin_genome_2 = origin_genome_long_2[2:-2]
                        mapped_strand_2 = "+"

                    elif nn==2: 							# RC-FW mapped to + strand:

#                        original_BS_1 = original_bs_reads_2[header]
#                        original_BS_2 = reverse_compl_seq(original_bs_reads_1[header])
                        FR = "+RF"
#                        mapped_location_1 += 1
#                        origin_genome_long_1 = my_gseq[mapped_location_1 - 2 - 1 : mapped_location_1 + g_len_1 + 2 - 1]
#                        origin_genome_1 = origin_genome_long_1[2:-2]
                        mapped_strand_1 = "+"

#                        mapped_location_2 += 1
#                        origin_genome_long_2 = my_gseq[mapped_location_2 - 2 - 1 : mapped_location_2 + g_len_2 + 2 - 1]
#                        origin_genome_2 = origin_genome_long_2[2:-2]
                        mapped_strand_2 = "+"


                    elif nn==3: 						# FW-RC mapped to - strand:
#                        original_BS_1=original_bs_reads_1[header]
#                        original_BS_2=reverse_compl_seq(original_bs_reads_2[header])

                        FR = "-FR"
                        mapped_location_1 = chr_length - mapped_location_1 - g_len_1
#                        origin_genome_long_1 = my_gseq[mapped_location_1 - 2 - 1 : mapped_location_1 + g_len_1 + 2 - 1]
#                        origin_genome_long_1 = reverse_compl_seq(origin_genome_long_1)
#                        origin_genome_1 = origin_genome_long_1[2:-2]
                        mapped_strand_1 = "-"

                        mapped_location_2 = chr_length - mapped_location_2 - g_len_2
#                        origin_genome_long_2 = reverse_compl_seq(my_gseq[mapped_location_2 - 2 - 1 : mapped_location_2 + g_len_2 + 2 - 1 ])
#                        origin_genome_2 = origin_genome_long_2[2:-2]
                        mapped_strand_2 = "-"

                    elif nn==4: 						# RC-FW mapped to - strand:
#                        original_BS_1=original_bs_reads_2[header]
#                        original_BS_2=reverse_compl_seq(original_bs_reads_1[header])

                        FR = "-RF"
                        mapped_location_1 = chr_length - mapped_location_1 - g_len_1
#                        origin_genome_long_1 = my_gseq[mapped_location_1 - 2 - 1 : mapped_location_1 + g_len_1 + 2 - 1]
#                        origin_genome_long_1 = reverse_compl_seq(origin_genome_long_1)
#                        origin_genome_1 = origin_genome_long_1[2:-2]
                        mapped_strand_1 = "-"

                        mapped_location_2 = chr_length - mapped_location_2 - g_len_2
#                        origin_genome_long_2 = reverse_compl_seq(my_gseq[mapped_location_2 - 2 - 1 : mapped_location_2 + g_len_2 + 2 - 1])
#                        origin_genome_2 = origin_genome_long_2[2:-2]
                        mapped_strand_2 = "-"




                    origin_genome_1, next_1, output_genome_1 = get_genomic_sequence(my_gseq, mapped_location_1, mapped_location_1 + g_len_1, mapped_strand_1)
                    origin_genome_2, next_2, output_genome_2 = get_genomic_sequence(my_gseq, mapped_location_2, mapped_location_2 + g_len_2, mapped_strand_2)

                    r_aln_1, g_aln_1 = cigar_to_alignment(cigar1, original_BS_1, origin_genome_1)
                    r_aln_2, g_aln_2 = cigar_to_alignment(cigar2, original_BS_2, origin_genome_2)


                    N_mismatch_1 = N_MIS(r_aln_1, g_aln_1) #+ original_BS_length_1 - (r_end_1 - r_start_1) # mismatches in the alignment + soft clipped nucleotides
                    N_mismatch_2 = N_MIS(r_aln_2, g_aln_2) #+ original_BS_length_2 - (r_end_2 - r_start_2) # mismatches in the alignment + soft clipped nucleotides


                    if max(N_mismatch_1, N_mismatch_2) <= int(max_mismatch_no) :
                        all_mapped_passed += 1
                        numbers_mapped_lst[nn-1] += 1
                        #---- unmapped -------------------------
                        del original_bs_reads_1[header]
                        del original_bs_reads_2[header]
                        #---------------------------------------

#                        output_genome_1 = origin_genome_long_1[0:2] + "_" + origin_genome_1 + "_" + origin_genome_long_1[-2:]
#                        output_genome_2 = origin_genome_long_2[0:2] + "_" + origin_genome_2 + "_" + origin_genome_long_2[-2:]

                        methy_1=methy_seq(r_aln_1, g_aln_1 + next_1)
                        methy_2=methy_seq(r_aln_2, g_aln_2 + next_2)

                        mC_lst, uC_lst = mcounts(methy_1, mC_lst, uC_lst)
                        mC_lst, uC_lst = mcounts(methy_2, mC_lst, uC_lst)


                        # 
                        #---XS FILTER----------------
                        #XS = 1 if "ZZZ" in methy.replace('-', '') else 0
                        XS_1 = 0
                        nCH_1 = methy_1.count('y') + methy_1.count('z')
                        nmCH_1 = methy_1.count('Y') + methy_1.count('Z')
                        if( (nmCH_1>XS_count) and nmCH_1/float(nCH_1+nmCH_1)>XS_pct ) :
                            XS_1 = 1
                        #XS = 1 if "ZZZ" in methy.replace('-', '') else 0
                        XS_2 = 0
                        nCH_2 = methy_2.count('y') + methy_2.count('z')
                        nmCH_2 = methy_2.count('Y') + methy_2.count('Z')
                        if( (nmCH_2>XS_count) and nmCH_2/float(nCH_2+nmCH_2)>XS_pct ) :
                            XS_2 = 1


                        outfile.store(header, N_mismatch_1, FR, mapped_chr, mapped_strand_1, mapped_location_1, cigar1, original_BS_1, methy_1, XS_1,
                            output_genome = output_genome_1, rnext = mapped_chr, pnext = mapped_location_2)
                        outfile.store(header, N_mismatch_2, FR, mapped_chr, mapped_strand_2, mapped_location_2, cigar2, original_BS_2, methy_2, XS_2,
                            output_genome = output_genome_2, rnext = mapped_chr, pnext = mapped_location_1)

            print "--> %s %s (%d/%d) " % (read_file_1, read_file_2, no_my_files, len(my_files))
            #----------------------------------------------------------------
            #	output unmapped pairs
            #----------------------------------------------------------------

            unmapped_lst=original_bs_reads_1.keys()
            unmapped_lst.sort()

#            for u in unmapped_lst:
#                outf_u1.write("%s\n"%original_bs_reads_1[u])
#                outf_u2.write("%s\n"%original_bs_reads_2[u])

            all_unmapped += len(unmapped_lst)


        if asktag=="N":

            #----------------------------------------------------------------
            outfile_1FCT= tmp_d('Trimed_FCT_1.fa'+random_id)
            outfile_2FCT= tmp_d('Trimed_FCT_2.fa'+random_id)

            try :
                read_inf=open(tmp_d(read_file_1),"r")
            except IOError :
                print "[Error] Cannot open file : %s", tmp_d(read_file_1)
                exit(-1)

            oneline=read_inf.readline()
            l=oneline.split()
            input_format=""

            #if len(l)==5: # old solexa format
            #	input_format="old Solexa Seq file"

            if oneline[0]=="@":	# Illumina GAII FastQ (Lister et al Nature 2009)
                input_format="FastQ"
                n_fastq=0
            elif len(l)==1 and oneline[0]!=">": 	# pure sequences
                input_format="_list of sequences"
            elif len(l)==11:	# Illumina GAII qseq file
                input_format="Illumina GAII qseq file"
            elif oneline[0]==">":	# fasta
                input_format="fasta"
                n_fasta=0

            read_inf.close()

            print "Detected data format: %s" % input_format

            #----------------------------------------------------------------
            read_file_list=[read_file_1,read_file_2]
            outfile_FCT_list=[outfile_1FCT,outfile_2FCT]
            n_list=[0,0]

            for f in range(2):
                read_file=read_file_list[f]
                outf_FCT=open(outfile_FCT_list[f],'w')
                original_bs_reads = original_bs_reads_lst[f]
                n=n_list[f]

                seq_id=""
                seq=""
                seq_ready="N"
                for line in fileinput.input(tmp_d(read_file)):
                    l=line.split()
                    if input_format=="old Solexa Seq file":
                        n+=1
                        seq_id=str(n)
                        seq_id=seq_id.zfill(12)
                        seq=l[4]
                        seq_ready="Y"
                    elif input_format=="_list of sequences":
                        n+=1
                        seq_id=str(n)
                        seq_id=seq_id.zfill(12)
                        seq=l[0]
                        seq_ready="Y"
                    elif input_format=="FastQ":
                        m_fastq=math.fmod(n_fastq,4)
                        n_fastq+=1
                        seq_ready="N"
                        if m_fastq==0:
                            n+=1
                            seq_id=str(n)
                            seq_id=seq_id.zfill(12)
                            seq=""
                        elif m_fastq==1:
                            seq=l[0]
                            seq_ready="Y"
                        else:
                            seq=""
                    elif input_format=="Illumina GAII qseq file":
                        n+=1
                        seq_id=str(n)
                        seq_id=seq_id.zfill(12)
                        seq=l[8]
                        seq_ready="Y"
                    elif input_format=="fasta":
                        m_fasta=math.fmod(n_fasta,2)
                        n_fasta+=1
                        seq_ready="N"
                        if m_fasta==0:
                            n+=1
                            seq_id=l[0][1:]
                            seq_id=seq_id.zfill(17)
                            seq=""
                        elif m_fasta==1:
                            seq=l[0]
                            seq_ready="Y"
                        else:
                            seq=""
                    #----------------------------------------------------------------
                    if seq_ready=="Y":
                        seq=seq[cut1-1:cut2] #<----------------------selecting 0..52 from 1..72  -e 52
                        seq=seq.upper()
                        seq=seq.replace(".","N")

                        #--striping BS adapter from 3' read --------------------------------------------------------------
                        if (adapterA !="") and (adapterB !=""):
                            signature=adapterA[:6]
                            if signature in seq:
                                signature_pos=seq.index(signature)
                                if seq[signature_pos:] in adapterA:
                                    seq=seq[:signature_pos]#+"".join(["N" for x in range(len(seq)-len(signature_pos))])
                                    all_trimmed+=1
                            else:
                                signature=adapterB[:6]
                                if signature in seq:
                                    #print seq_id,seq,signature;
                                    signature_pos=seq.index(signature)
                                    if seq[signature_pos:] in adapterB:
                                        seq=seq[:signature_pos]#+"".join(["N" for x in range(len(seq)-len(signature_pos))])
                                        all_trimmed+=1

                        if len(seq) <= 4:
                            seq = "N" * (cut2-cut1+1)
                        #---------  trimmed_raw_BS_read  ------------------
                        original_bs_reads[seq_id] = seq

                        #---------  FW_C2T  ------------------
                        if f==0:
                            outf_FCT.write('>%s\n%s\n'% (seq_id, seq.replace("C","T")))
                        elif f==1:
                            outf_FCT.write('>%s\n%s\n'% (seq_id, reverse_compl_seq(seq).replace("C","T")))


                n_list[f]=n
                outf_FCT.close()

                fileinput.close()


            #print "All input end 1: %d , end 2: %d "%(n_list[0],n_list[1]);
            all_raw_reads+=n

            #--------------------------------------------------------------------------------
            # Bowtie mapping
            #--------------------------------------------------------------------------------
            WC2T_fr=tmp_d("W_C2T_fr_m"+max_mismatch_no+".mapping"+random_id)
            CC2T_fr=tmp_d("C_C2T_fr_m"+max_mismatch_no+".mapping"+random_id)

            run_in_parallel([ aligner_command % {'reference_genome' : os.path.join(db_path,'W_C2T'),
                                         'input_file_1' : outfile_1FCT,
                                         'input_file_2' : outfile_2FCT,
                                         'output_file' : WC2T_fr},

                              aligner_command % {'reference_genome' : os.path.join(db_path,'C_C2T'),
                                         'input_file_1' : outfile_1FCT,
                                         'input_file_2' : outfile_2FCT,
                                         'output_file' : CC2T_fr} ])


            delete_files(outfile_1FCT, outfile_2FCT)

            #--------------------------------------------------------------------------------
            # Post processing
            #--------------------------------------------------------------------------------

            FW_C2T_fr_U, FW_C2T_fr_R = extract_mapping(WC2T_fr)
            RC_C2T_fr_U, RC_C2T_fr_R = extract_mapping(CC2T_fr)

            #----------------------------------------------------------------
            # get uniq-hit reads
            #----------------------------------------------------------------
            Union_set = set(FW_C2T_fr_U.iterkeys()) | set(RC_C2T_fr_U.iterkeys())

            Unique_FW_fr_C2T = set() # +
            Unique_RC_fr_C2T = set() # -
            Multiple_hits=set()


            for x in Union_set:
                _list = []
                for d in [FW_C2T_fr_U, RC_C2T_fr_U]:
                    mis_lst = d.get(x,[99])
                    mis = int(mis_lst[0])
                    _list.append(mis)
                for d in [FW_C2T_fr_R, RC_C2T_fr_R]:
                    mis = d.get(x,99)
                    _list.append(mis)
                mini = min(_list)
                if _list.count(mini) == 1:
                    mini_index = _list.index(mini)
                    if mini_index == 0:
                        Unique_FW_fr_C2T.add(x)
                    elif mini_index == 1:
                        Unique_RC_fr_C2T.add(x)
                    else :
                        Multiple_hits.add(x)
                else :
                    Multiple_hits.add(x)

           # write reads rejected by Multiple Hits to file
            if show_multiple_hit :
                outf_MH=open("Multiple_hit.fa",'w')
                for i in Multiple_hits :
                    outf_MH.write(">%s\n" % i)
                    outf_MH.write("%s\n" % original_bs_reads[i])
                outf_MH.close()

            FW_C2T_fr_uniq_lst=[[FW_C2T_fr_U[u][1],u] for u in Unique_FW_fr_C2T]
            RC_C2T_fr_uniq_lst=[[RC_C2T_fr_U[u][1],u] for u in Unique_RC_fr_C2T]

            FW_C2T_fr_uniq_lst.sort()
            RC_C2T_fr_uniq_lst.sort()
            FW_C2T_fr_uniq_lst=[x[1] for x in FW_C2T_fr_uniq_lst]
            RC_C2T_fr_uniq_lst=[x[1] for x in RC_C2T_fr_uniq_lst]
            #----------------------------------------------------------------

            numbers_premapped_lst[0]+=len(Unique_FW_fr_C2T)
            numbers_premapped_lst[1]+=len(Unique_RC_fr_C2T)


            #----------------------------------------------------------------

            nn = 0
            for ali_unique_lst, ali_dic in [(FW_C2T_fr_uniq_lst,FW_C2T_fr_U), (RC_C2T_fr_uniq_lst,RC_C2T_fr_U)]:
                nn += 1
                mapped_chr0 = ""
                for header in ali_unique_lst:
                    _, mapped_chr, mapped_location_1, cigar1, mapped_location_2, cigar2 = ali_dic[header]


                    #-------------------------------------
                    if mapped_chr != mapped_chr0:
                        my_gseq = deserialize(db_d(mapped_chr))
                        chr_length = len(my_gseq)
                        mapped_chr0 = mapped_chr
                    #-------------------------------------

                    original_BS_1 = original_bs_reads_1[header]
                    original_BS_2 = reverse_compl_seq(original_bs_reads_2[header])

                    r_start_1, r_end_1, g_len_1 = get_read_start_end_and_genome_length(cigar1)
                    r_start_2, r_end_2, g_len_2 = get_read_start_end_and_genome_length(cigar2)

                    all_mapped += 1

                    if nn == 1: 							# FW-RC mapped to + strand:
                        FR = "+FR"
#                        mapped_location_1 += 1
#                        origin_genome_long_1 = my_gseq[mapped_location_1 - 2 - 1 : mapped_location_1 + g_len_1 + 2 - 1]
#                        origin_genome_1 = origin_genome_long_1[2:-2]
                        mapped_strand_1 = "+"

#                        mapped_location_2 += 1
#                        origin_genome_long_2 = my_gseq[mapped_location_2 - 2 - 1 : mapped_location_2 + g_len_2 + 2 - 1]
#                        origin_genome_2 = origin_genome_long_2[2:-2]
                        mapped_strand_2 = "+"

                    elif nn == 2: 						# FW-RC mapped to - strand:

                        FR="-FR"
                        mapped_location_1 = chr_length - mapped_location_1 - g_len_1
#                        origin_genome_long_1 = my_gseq[mapped_location_1 - 2 - 1 : mapped_location_1 + g_len_1 + 2 - 1]
#                        origin_genome_long_1 = reverse_compl_seq(origin_genome_long_1)
#                        origin_genome_1 = origin_genome_long_1[2:-2]
                        mapped_strand_1 = "-"

                        mapped_location_2 = chr_length - mapped_location_2 - g_len_2
#                        origin_genome_long_2 = reverse_compl_seq(my_gseq[mapped_location_2 - 2 - 1 : mapped_location_2 + g_len_2 + 2 - 1])
#                        origin_genome_2 = origin_genome_long_2[2:-2]
                        mapped_strand_2 = "-"



                    origin_genome_1, next_1, output_genome_1 = get_genomic_sequence(my_gseq, mapped_location_1, mapped_location_1 + g_len_1, mapped_strand_1)
                    origin_genome_2, next_2, output_genome_2 = get_genomic_sequence(my_gseq, mapped_location_2, mapped_location_2 + g_len_2, mapped_strand_2)


                    r_aln_1, g_aln_1 = cigar_to_alignment(cigar1, original_BS_1, origin_genome_1)
                    r_aln_2, g_aln_2 = cigar_to_alignment(cigar2, original_BS_2, origin_genome_2)

                    N_mismatch_1 = N_MIS(r_aln_1, g_aln_1) #+ original_BS_length_1 - (r_end_1 - r_start_1) # mismatches in the alignment + soft clipped nucleotides
                    N_mismatch_2 = N_MIS(r_aln_2, g_aln_2) #+ original_BS_length_2 - (r_end_2 - r_start_2) # mismatches in the alignment + soft clipped nucleotides

                    if max(N_mismatch_1, N_mismatch_2) <= int(max_mismatch_no):

                        numbers_mapped_lst[nn-1] += 1
                        all_mapped_passed += 1

                        #---- unmapped -------------------------
                        del original_bs_reads_1[header]
                        del original_bs_reads_2[header]

                        #---------------------------------------
#                        output_genome_1 = origin_genome_long_1[0:2] + "_" + origin_genome_1 + "_" + origin_genome_long_1[-2:]
#                        output_genome_2 = origin_genome_long_2[0:2] + "_" + origin_genome_2 + "_" + origin_genome_long_2[-2:]
                        methy_1=methy_seq(r_aln_1, g_aln_1 + next_1)
                        methy_2=methy_seq(r_aln_2, g_aln_2 + next_2)
                        mC_lst,uC_lst = mcounts(methy_1, mC_lst, uC_lst)
                        mC_lst,uC_lst = mcounts(methy_2, mC_lst, uC_lst)

                        # 
                        #---XS FILTER----------------
                        #XS = 1 if "ZZZ" in methy.replace('-', '') else 0
                        XS_1 = 0
                        nCH_1 = methy_1.count('y') + methy_1.count('z')
                        nmCH_1 = methy_1.count('Y') + methy_1.count('Z')
                        if( (nmCH_1>XS_count) and nmCH_1/float(nCH_1+nmCH_1)>XS_pct ) :
                            XS_1 = 1
                        #XS = 1 if "ZZZ" in methy.replace('-', '') else 0
                        XS_2 = 0
                        nCH_2 = methy_2.count('y') + methy_2.count('z')
                        nmCH_2 = methy_2.count('Y') + methy_2.count('Z')
                        if( (nmCH_2>XS_count) and nmCH_2/float(nCH_2+nmCH_2)>XS_pct ) :
                            XS_2 = 1


                        outfile.store(header, N_mismatch_1, FR, mapped_chr, mapped_strand_1, mapped_location_1, cigar1, original_BS_1, methy_1, XS_1,
                                      output_genome = output_genome_1, rnext = mapped_chr, pnext = mapped_location_2)
                        outfile.store(header, N_mismatch_2, FR, mapped_chr, mapped_strand_2, mapped_location_2, cigar2, original_BS_2, methy_2, XS_2,
                                      output_genome = output_genome_2, rnext = mapped_chr, pnext = mapped_location_1)

            print "--> %s %s (%d/%d) " % (read_file_1, read_file_2, no_my_files, len(my_files))
            #----------------------------------------------------------------
            #	output unmapped pairs
            #----------------------------------------------------------------

            unmapped_lst=original_bs_reads_1.keys()
            unmapped_lst.sort()

#            for u in unmapped_lst:
#                outf_u1.write("%s\n"%(original_bs_reads_1[u]))
#                outf_u2.write("%s\n"%(original_bs_reads_2[u]) )

            all_unmapped+=len(unmapped_lst)


    #==================================================================================================
#    outf.close()
#
#    outf_u1.close()
#    outf_u2.close()
    delete_files(tmp_path)

    logm("-------------------------------- " )
    logm("O Number of raw BS-read pairs: %d ( %d bp)"%(all_raw_reads,cut2-cut1+1) )
    logm("O Number of ends trimmed for adapter: %d"% all_trimmed+"\n")

    if all_raw_reads >0:

        logm("O Number of unique-hits read pairs for post-filtering: %d" % all_mapped + "\n")
        if asktag=="Y":
            logm("O -- %7d FW-RC pairs mapped to Watson strand (before post-filtering)"%(numbers_premapped_lst[0]) )
            logm("O -- %7d RC-FW pairs mapped to Watson strand (before post-filtering)"%(numbers_premapped_lst[1]) )
            logm("O -- %7d FW-RC pairs mapped to Crick strand (before post-filtering)"%(numbers_premapped_lst[2]) )
            logm("O -- %7d RC-FW pairs mapped to Crick strand (before post-filtering)"%(numbers_premapped_lst[3]) )
        elif asktag=="N":
            logm("O -- %7d FW-RC pairs mapped to Watson strand (before post-filtering)"%(numbers_premapped_lst[0]) )
            logm("O -- %7d FW-RC pairs mapped to Crick strand (before post-filtering)"%(numbers_premapped_lst[1]) )


        logm("O --- Number of reads rejected because of multiple hits: %d\n" % len(Multiple_hits) )
        logm("O --- %d uniquely aligned pairs, where each end has mismatches <= %s"%(all_mapped_passed, max_mismatch_no) )
        if asktag=="Y":
            logm("O ----- %7d FW-RC pairs mapped to Watson strand"%(numbers_mapped_lst[0]) )
            logm("O ----- %7d RC-FW pairs mapped to Watson strand"%(numbers_mapped_lst[1]) )
            logm("O ----- %7d FW-RC pairs mapped to Crick strand"%(numbers_mapped_lst[2]) )
            logm("O ----- %7d RC-FW pairs mapped to Crick strand"%(numbers_mapped_lst[3]) )
        elif asktag=="N":
            logm("O ----- %7d FW-RC pairs mapped to Watson strand"%(numbers_mapped_lst[0]) )
            logm("O ----- %7d FW-RC pairs mapped to Crick strand"%(numbers_mapped_lst[1]) )

        logm("O Mappability= %1.4f%%"%(100*float(all_mapped_passed)/all_raw_reads) )
        logm("O Unmapped read pairs: %d"% all_unmapped+"\n")


        n_CG=mC_lst[0]+uC_lst[0]
        n_CHG=mC_lst[1]+uC_lst[1]
        n_CHH=mC_lst[2]+uC_lst[2]

        logm("-------------------------------- " )
        logm("Methylated C in mapped reads " )
        logm(" mCG %1.3f%%"%((100*float(mC_lst[0])/n_CG) if n_CG != 0 else 0) )
        logm(" mCHG %1.3f%%"%((100*float(mC_lst[1])/n_CHG) if n_CHG != 0 else 0) )
        logm(" mCHH %1.3f%%"%((100*float(mC_lst[2])/n_CHH) if n_CHH != 0 else 0) )

    logm("----------------------------------------------" )
    logm("------------------- END ----------------------" )
    elapsed("=== END %s %s ===" % (main_read_file_1, main_read_file_2))
