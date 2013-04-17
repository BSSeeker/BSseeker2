import fileinput, random, math, os.path
from bs_index.rrbs_build import FWD_MAPPABLE_REGIONS, REV_MAPPABLE_REGIONS
from bs_utils.utils import *

from bs_align.bs_single_end import extract_mapping
from bs_align_utils import *

def my_mappable_region(chr_regions, mapped_location, FR): # start_position (first C), end_position (last G), serial, sequence
    #print len(chr_regions)
    out_serial=0
    out_start=-1
    out_end=-1
    if FR == "+FW":
        my_location=str(mapped_location-2)
        if my_location in chr_regions:
            my_lst = chr_regions[my_location]
            out_start = int(my_location)
            out_end = my_lst[0]
            out_serial = my_lst[1]
        #else :
        #    print "[For debug]: +FW location %s cannot be found" % my_location

    elif FR == "-FW":
        my_location = str(mapped_location-1)
        if my_location in chr_regions:
            my_lst = chr_regions[my_location]
            out_end = int(my_location)
            out_start = my_lst[0]
            out_serial = my_lst[1]
        #else :
        #    print "[For debug]: -FW location %s cannot be found" % my_location

    return out_serial, out_start, out_end



#----------------------------------------------------------------

def bs_rrbs(main_read_file, asktag, mytag, adapter_file, cut_s, cut_e, no_small_lines, max_mismatch_no,
            aligner_command, db_path, tmp_path, outfile, XS_pct, XS_count, adapter_mismatch):
    #----------------------------------------------------------------

    mytag_lst = mytag.split("/")
    #----------------------------------------------------------------

    # helper method to join fname with tmp_path
    tmp_d = lambda fname: os.path.join(tmp_path, fname)
    db_d = lambda fname:  os.path.join(db_path, fname)

    #----------------------------------------------------------------
    adapter_seq=""
    if adapter_file:
        try :
            adapter_inf = open(adapter_file,"r")
            adapter_seq = adapter_inf.readline().strip()
            adapter_inf.close()
        except IOError:
            print "[Error] Cannot find adapter file : %s !" % adapter_file
            exit(-1)

    logm("I Read filename: %s" % main_read_file)
    logm("I Starting Msp-1 tag: %s" % mytag )
    logm("I The last cycle (for mapping): %d" % cut_e )
    logm("I Bowtie path: %s" % aligner_command )
    logm("I Reference genome library path: %s" % db_path )
    logm("I Number of mismatches allowed: %s" % max_mismatch_no)
    logm("I Adapter seq: %s" % adapter_seq)
    logm("----------------------------------------------")

    #----------------------------------------------------------------
    all_raw_reads=0
    all_tagged=0
    all_tagged_trimed=0
    all_mapped=0
    all_mapped_passed=0
    n_mytag_lst={}
    for x in mytag_lst:
        n_mytag_lst[x]=0

    mC_lst=[0,0,0]
    uC_lst=[0,0,0]

    no_my_files=0

    num_mapped_FW_C2T = 0
    num_mapped_RC_C2T = 0
    num_mapped_FW_G2A = 0
    num_mapped_RC_G2A = 0

    #===============================================
    # directional sequencing
    #===============================================

    if asktag=="N" :
        #----------------------------------------------------------------
        logm("== Start mapping ==")

        input_fname = os.path.split(main_read_file)[1]
        for read_file in isplit_file(main_read_file, tmp_d(input_fname)+'-s-', no_small_lines):

            logm("Processing read file: %s" % read_file)
            original_bs_reads = {}
            no_my_files+=1
            random_id = ".tmp-"+str(random.randint(1000000,9999999))
            outfile2=tmp_d('Trimed_C2T.fa'+random_id)

            outf2=open(outfile2,'w')

            #--- Checking input format ------------------------------------------
            try :
                read_inf=open(read_file,"r")
            except IOError:
                print "[Error] Cannot open input file : %s" % read_file
                exit(-1)

            oneline=read_inf.readline()
            l=oneline.split()
            n_fastq=0
            n_fasta=0
            input_format=""
            if oneline[0]=="@":	# FastQ
                input_format="fastq"
            #    n_fastq=0
            elif len(l)==1 and oneline[0]!=">": # pure sequences
                input_format="seq"
            elif len(l)==11: # Illumina qseq
                input_format="qseq"
            elif oneline[0]==">": # fasta
                input_format="fasta"
            #    n_fasta=0
            read_inf.close()

            #----------------------------------------------------------------
            seq_id=""
            seq=""
            seq_ready=0
            for line in fileinput.input(read_file):
                l=line.split()

                if input_format=="seq":
                    all_raw_reads+=1
                    seq_id=str(all_raw_reads)
                    seq_id=seq_id.zfill(12)
                    seq=l[0]
                    seq_ready="Y"
                elif input_format=="fastq":
                    m_fastq=math.fmod(n_fastq,4)
                    n_fastq+=1
                    seq_ready="N"
                    if m_fastq==0:
                        all_raw_reads+=1
                        seq_id=str(all_raw_reads)
                        seq_id=seq_id.zfill(12)
                        seq=""
                    elif m_fastq==1:
                        seq=l[0]
                        seq_ready="Y"
                    else:
                        seq=""
                elif input_format=="qseq":
                    all_raw_reads+=1
                    seq_id=str(all_raw_reads)
                    seq_id=seq_id.zfill(12)
                    seq=l[8]
                    seq_ready="Y"
                elif input_format=="fasta":
                    m_fasta=math.fmod(n_fasta,2)
                    n_fasta+=1
                    seq_ready="N"
                    if m_fasta==0:
                        all_raw_reads+=1
                        seq_id=l[0][1:]
                        seq=""
                    elif m_fasta==1:
                        seq=l[0]
                        seq_ready="Y"
                    else:
                        seq=""
                #---------------------------------------------------------------
                if seq_ready=="Y":
                    # Normalize the characters
                    seq=seq.upper().replace(".","N")

                    if seq[0:3] in mytag_lst :
                        all_tagged+=1
                        n_mytag_lst[seq[0:3]]+=1

                    seq = seq[(cut_s-1):cut_e] # cut_s start from 1 cycle by default

                    #-- Trimming adapter sequence ---
                    if adapter_seq != "" :
                        new_read = RemoveAdapter(seq, adapter_seq, adapter_mismatch)
                        if len(new_read) < len(seq) :
                            all_tagged_trimed += 1
                        seq = new_read
                    if len(seq) <= 4 :
                        seq = "N" * (cut_e - cut_s)

                    # all reads will be considered, regardless of tags
                    #---------  trimmed_raw_BS_read and qscore ------------------
                    original_bs_reads[seq_id] = seq
                    #---------  FW_C2T  ------------------
                    outf2.write('>%s\n%s\n'%(seq_id, seq.replace('C', 'T')))
            fileinput.close()

            outf2.close()

            delete_files(read_file)
            logm("Processing input is done")
            #--------------------------------------------------------------------------------

            # mapping
            #--------------------------------------------------------------------------------
            WC2T=tmp_d("W_C2T_m"+max_mismatch_no+".mapping"+random_id)
            CC2T=tmp_d("C_C2T_m"+max_mismatch_no+".mapping"+random_id)

            run_in_parallel([ aligner_command % {'reference_genome' : os.path.join(db_path,'W_C2T'),
                                                   'input_file' : outfile2,
                                                   'output_file' : WC2T},
                              aligner_command % {'reference_genome' : os.path.join(db_path,'C_C2T'),
                                                   'input_file' : outfile2,
                                                   'output_file' : CC2T} ])

            logm("Aligning reads is done")

            delete_files(outfile2)

            #--------------------------------------------------------------------------------
            # Post processing
            #--------------------------------------------------------------------------------

            FW_C2T_U,FW_C2T_R=extract_mapping(WC2T)
            RC_C2T_U,RC_C2T_R=extract_mapping(CC2T)
            logm("Extracting alignments is done")

            #----------------------------------------------------------------
            # get uniq-hit reads
            #----------------------------------------------------------------
            Union_set=set(FW_C2T_U.iterkeys()) | set(RC_C2T_U.iterkeys())

            Unique_FW_C2T=set() # +
            Unique_RC_C2T=set() # -


            for x in Union_set:
                _list=[]
                for dx in [FW_C2T_U, RC_C2T_U]:
                    mis_lst=dx.get(x,[99])
                    mis=int(mis_lst[0])
                    _list.append(mis)
                for dx in [FW_C2T_R, RC_C2T_R]:
                    mis=dx.get(x,99)
                    _list.append(mis)
                mini=min(_list)
                if _list.count(mini)==1:
                    mini_index=_list.index(mini)
                    if mini_index==0:
                        Unique_FW_C2T.add(x)
                    elif mini_index==1:
                        Unique_RC_C2T.add(x)

            del Union_set
            del FW_C2T_R
            del RC_C2T_R

            FW_uniq_lst=[[FW_C2T_U[u][1],u] for u in Unique_FW_C2T]
            RC_uniq_lst=[[RC_C2T_U[u][1],u] for u in Unique_RC_C2T]
            FW_uniq_lst.sort()
            RC_uniq_lst.sort()
            FW_uniq_lst=[x[1] for x in FW_uniq_lst]
            RC_uniq_lst=[x[1] for x in RC_uniq_lst]

            del Unique_FW_C2T
            del Unique_RC_C2T

            #----------------------------------------------------------------
            # Post-filtering reads
            # ---- FW  ----
            FW_regions = dict()
            gseq = dict()
            chr_length = dict()
            for header in FW_uniq_lst :
                _, mapped_chr, mapped_location, cigar = FW_C2T_U[header]
                original_BS = original_bs_reads[header]
                if mapped_chr not in FW_regions :
                    FW_regions[mapped_chr] = deserialize(db_d(FWD_MAPPABLE_REGIONS(mapped_chr)))
                if mapped_chr not in gseq :
                    gseq[mapped_chr] = deserialize(db_d(mapped_chr))
                    chr_length[mapped_chr] = len(gseq[mapped_chr])

                r_start, r_end, g_len = get_read_start_end_and_genome_length(cigar)
                all_mapped+=1
                FR = "+FW"
                mapped_strand = "+"
                origin_genome, next2bp, output_genome = get_genomic_sequence(gseq[mapped_chr],
                                                                             mapped_location,
                                                                             mapped_location + g_len,
                                                                             mapped_strand)
                checking_first_C = (output_genome[1:6] == "C_CGG")
                r_aln, g_aln = cigar_to_alignment(cigar, original_BS, origin_genome)

                if len(r_aln) == len(g_aln) :
                    my_region_serail = 0
                    if checking_first_C :
                        my_region_serial, my_region_start, my_region_end = my_mappable_region(FW_regions[mapped_chr],
                                                                                         mapped_location + 1,
                                                                                         FR)
                    if my_region_serial == 0 : # still be 0
                        # for some cases, read has no tags; searching the upstream sequence for CCGG tags
                        # print "[For debug]: FW read has no tags"
                        upstream, _, _ = get_genomic_sequence(gseq[mapped_chr], mapped_location - 500,
                                                              mapped_location + 4, '-')
                        # print "[For debug]: upstream=", upstream
                        # print "[For debug]: mapped_locaiton: %d, CCGG_pos: %d" %(mapped_location, upstream.find('CCGG'))
                        my_region_serial, my_region_start, my_region_end = my_mappable_region(FW_regions[mapped_chr],
                                                              mapped_location + 2 - upstream.find('CCGG'), FR)
                        #if my_region_serial == 0 :
                        #    print "[For debug]: chr=", mapped_chr
                        #    print "[For debug]: +FW read still can not find fragment serial"
                        # Tip: sometimes "my_region_serial" is still 0 ...


                    N_mismatch = N_MIS(r_aln, g_aln)
                    if N_mismatch <= int(max_mismatch_no) :
                        all_mapped_passed += 1
                        methy = methy_seq(r_aln, g_aln + next2bp)
                        mC_lst, uC_lst = mcounts(methy, mC_lst, uC_lst)
                        #---XS FILTER----------------
                        XS = 0
                        nCH = methy.count('y') + methy.count('z')
                        nmCH = methy.count('Y') + methy.count('Z')
                        if( (nmCH>XS_count) and nmCH/float(nCH+nmCH)>XS_pct ) :
                            XS = 1
                        num_mapped_FW_C2T += 1
                        outfile.store(header, N_mismatch, FR, mapped_chr, mapped_strand,
                                      mapped_location, cigar, original_BS, methy, XS,
                                      output_genome = output_genome,
                                      rrbs = True,
                                      my_region_serial = my_region_serial,
                                      my_region_start = my_region_start,
                                      my_region_end = my_region_end)
                else :
                    print "[For debug]: reads not in same lengths"


            # ---- RC ----
            RC_regions = dict()
            for header in RC_uniq_lst :
                _, mapped_chr, mapped_location, cigar = RC_C2T_U[header]
                original_BS = original_bs_reads[header]
                if mapped_chr not in RC_regions :
                    RC_regions[mapped_chr] = deserialize(db_d(REV_MAPPABLE_REGIONS(mapped_chr)))
                if mapped_chr not in gseq :
                    gseq[mapped_chr] = deserialize(db_d(mapped_chr))
                    chr_length[mapped_chr] = len(gseq[mapped_chr])

                r_start, r_end, g_len = get_read_start_end_and_genome_length(cigar)
                mapped_location = chr_length[mapped_chr] - mapped_location - g_len
                all_mapped+=1
                FR = "-FW"
                mapped_strand = "-"
                origin_genome, next2bp, output_genome = get_genomic_sequence(gseq[mapped_chr],
                                                                             mapped_location,
                                                                             mapped_location + g_len,
                                                                             mapped_strand)
                checking_first_C = (output_genome[1:6] == "C_CGG")
                r_aln, g_aln = cigar_to_alignment(cigar, original_BS, origin_genome)

                if len(r_aln) == len(g_aln) : # and checking_first_C:
                    my_region_serial = 0
                    if checking_first_C :
                        my_region_serial, my_region_start, my_region_end = my_mappable_region(RC_regions[mapped_chr],
                                                                                         mapped_location + g_len + 1,
                                                                                         FR)
                    if my_region_serial == 0 : # still be 0
                        # for some cases, read has no tags; searching the upstream sequence for CCGG tags
                        # print "[For debug]: RC Read has no tags"
                        upstream, _, _ = get_genomic_sequence(gseq[mapped_chr], mapped_location - 3,
                                                              mapped_location + 500, '+')
                        #print "[For debug]: upstream = ", upstream
                        #print "[For debug]: mapped_location = %d, CCGG_pos: %d" % (mapped_location, upstream.find('CCGG'))
                        my_region_serial, my_region_start, my_region_end = my_mappable_region(RC_regions[mapped_chr],
                                        mapped_location + 1 + upstream.find('CCGG'), FR)
                        #if my_region_serial == 0 :
                        #    print "[For debug]: chr=", mapped_chr
                        #    print "[For debug]: -FW read still cannot find fragment serial"


                    N_mismatch = N_MIS(r_aln, g_aln)
                    if N_mismatch <= int(max_mismatch_no) :
                        all_mapped_passed += 1
                        methy = methy_seq(r_aln, g_aln + next2bp)
                        mC_lst, uC_lst = mcounts(methy, mC_lst, uC_lst)
                        #---XS FILTER----------------
                        XS = 0
                        nCH = methy.count('y') + methy.count('z')
                        nmCH = methy.count('Y') + methy.count('Z')
                        if( (nmCH>XS_count) and nmCH/float(nCH+nmCH)>XS_pct ) :
                            XS = 1
                        num_mapped_RC_C2T += 1
                        outfile.store(header, N_mismatch, FR, mapped_chr, mapped_strand,
                                      mapped_location, cigar, original_BS, methy, XS,
                                      output_genome = output_genome,
                                      rrbs = True,
                                      my_region_serial = my_region_serial,
                                      my_region_start = my_region_start,
                                      my_region_end = my_region_end)
                else :
                    print "[For debug]: reads not in same lengths"


            # Finished both FW and RC
            logm("Done: %s (%d) \n" % (read_file, no_my_files))
            print "--> %s (%d) "%(read_file, no_my_files)
            del original_bs_reads
            delete_files(WC2T, CC2T)

    # End of directional library


    # ====================================================
    #  directional library
    # ====================================================

    elif asktag=="Y" :
        #----------------------------------------------------------------
        logm("== Start mapping ==")

        input_fname = os.path.split(main_read_file)[1]
        for read_file in isplit_file(main_read_file, tmp_d(input_fname)+'-s-', no_small_lines):

            logm("Processing read file: %s" % read_file)
            original_bs_reads = {}
            no_my_files+=1
            random_id = ".tmp-"+str(random.randint(1000000,9999999))
            outfile2=tmp_d('Trimed_C2T.fa'+random_id)
            outfile3=tmp_d('Trimed_G2A.fa'+random_id)

            outf2=open(outfile2,'w')
            outf3=open(outfile3,'w')

            #--- Checking input format ------------------------------------------
            try :
                read_inf=open(read_file,"r")
            except IOError:
                print "[Error] Cannot open input file : %s" % read_file
                exit(-1)

            oneline=read_inf.readline()
            l=oneline.split()
            n_fastq=0
            n_fasta=0
            input_format=""
            if oneline[0]=="@":	# FastQ
                input_format="fastq"
            elif len(l)==1 and oneline[0]!=">": # pure sequences
                input_format="seq"
            elif len(l)==11: # Illumina qseq
                input_format="qseq"
            elif oneline[0]==">": # fasta
                input_format="fasta"
            read_inf.close()

            #----------------------------------------------------------------
            seq_id=""
            seq=""
            seq_ready=0
            for line in fileinput.input(read_file):
                l=line.split()

                if input_format=="seq":
                    all_raw_reads+=1
                    seq_id=str(all_raw_reads)
                    seq_id=seq_id.zfill(12)
                    seq=l[0]
                    seq_ready="Y"
                elif input_format=="fastq":
                    m_fastq=math.fmod(n_fastq,4)
                    n_fastq+=1
                    seq_ready="N"
                    if m_fastq==0:
                        all_raw_reads+=1
                        seq_id=str(all_raw_reads)
                        seq_id=seq_id.zfill(12)
                        seq=""
                    elif m_fastq==1:
                        seq=l[0]
                        seq_ready="Y"
                    else:
                        seq=""
                elif input_format=="qseq":
                    all_raw_reads+=1
                    seq_id=str(all_raw_reads)
                    seq_id=seq_id.zfill(12)
                    seq=l[8]
                    seq_ready="Y"
                elif input_format=="fasta":
                    m_fasta=math.fmod(n_fasta,2)
                    n_fasta+=1
                    seq_ready="N"
                    if m_fasta==0:
                        all_raw_reads+=1
                        seq_id=l[0][1:]
                        seq=""
                    elif m_fasta==1:
                        seq=l[0]
                        seq_ready="Y"
                    else:
                        seq=""
                    #---------------------------------------------------------------
                if seq_ready=="Y":
                    # Normalize the characters
                    seq=seq.upper().replace(".","N")

                    if seq[0:3] in mytag_lst :
                        all_tagged+=1
                        n_mytag_lst[seq[0:3]]+=1

                    seq = seq[(cut_s-1):cut_e] # cut_s start from 1 cycle by default

                    #-- Trimming adapter sequence ---
                    if adapter_seq != "" :
                        new_read = RemoveAdapter(seq, adapter_seq, adapter_mismatch)
                        if len(new_read) < len(seq) :
                            all_tagged_trimed += 1
                        seq = new_read
                    if len(seq) <= 4 :
                        seq = "N" * (cut_e - cut_s)

                    # all reads will be considered, regardless of tags
                    #---------  trimmed_raw_BS_read and qscore ------------------
                    original_bs_reads[seq_id] = seq
                    #---------  FW_C2T  ------------------
                    outf2.write('>%s\n%s\n'%(seq_id, seq.replace('C', 'T')))
                    #---------  RC_G2A  ------------------
                    outf3.write('>%s\n%s\n' % (seq_id, seq.replace("G","A")))
            fileinput.close()

            outf2.close()

            delete_files(read_file)
            logm("Processing input is done")
            #--------------------------------------------------------------------------------

            # mapping
            #--------------------------------------------------------------------------------
            WC2T=tmp_d("W_C2T_m"+max_mismatch_no+".mapping"+random_id)
            CC2T=tmp_d("C_C2T_m"+max_mismatch_no+".mapping"+random_id)
            WG2A=tmp_d("W_G2A_m"+max_mismatch_no+".mapping"+random_id)
            CG2A=tmp_d("C_G2A_m"+max_mismatch_no+".mapping"+random_id)

            run_in_parallel([ aligner_command % {'reference_genome' : os.path.join(db_path,'W_C2T'),
                                                 'input_file' : outfile2,
                                                 'output_file' : WC2T},
                              aligner_command % {'reference_genome' : os.path.join(db_path,'C_C2T'),
                                                 'input_file' : outfile2,
                                                 'output_file' : CC2T},
                              aligner_command % {'reference_genome' : os.path.join(db_path,'W_G2A'),
                                                 'input_file' : outfile3,
                                                 'output_file' : WG2A},
                              aligner_command % {'reference_genome' : os.path.join(db_path,'C_G2A'),
                                                 'input_file' : outfile3,
                                                 'output_file' : CG2A} ])

            logm("Aligning reads is done")

            delete_files(outfile2)

            #--------------------------------------------------------------------------------
            # Post processing
            #--------------------------------------------------------------------------------

            FW_C2T_U,FW_C2T_R=extract_mapping(WC2T)
            RC_G2A_U,RC_G2A_R=extract_mapping(CG2A)

            FW_G2A_U,FW_G2A_R=extract_mapping(WG2A)
            RC_C2T_U,RC_C2T_R=extract_mapping(CC2T)


            logm("Extracting alignments is done")

            #----------------------------------------------------------------
            # get uniq-hit reads
            #----------------------------------------------------------------
            Union_set=set(FW_C2T_U.iterkeys()) | set(RC_G2A_U.iterkeys()) | set(FW_G2A_U.iterkeys()) | set(RC_C2T_U.iterkeys())

            Unique_FW_C2T=set() # +
            Unique_RC_G2A=set() # +
            Unique_FW_G2A=set() # -
            Unique_RC_C2T=set() # -

            for x in Union_set:
                _list=[]
                for dx in [FW_C2T_U, RC_G2A_U, FW_G2A_U, RC_C2T_U]:
                    mis_lst=dx.get(x,[99])
                    mis=int(mis_lst[0])
                    _list.append(mis)
                for dx in [FW_C2T_R, RC_G2A_R, FW_G2A_R, RC_C2T_R]:
                    mis=dx.get(x,99)
                    _list.append(mis)
                mini=min(_list)
                if _list.count(mini) == 1:
                    mini_index=_list.index(mini)
                    if mini_index == 0:
                        Unique_FW_C2T.add(x)
                    elif mini_index == 1:
                        Unique_RC_G2A.add(x)
                    elif mini_index == 2:
                        Unique_FW_G2A.add(x)
                    elif mini_index == 3:
                        Unique_RC_C2T.add(x)



            del Union_set
            del FW_C2T_R
            del FW_G2A_R
            del RC_C2T_R
            del RC_G2A_R

            FW_C2T_uniq_lst=[[FW_C2T_U[u][1],u] for u in Unique_FW_C2T]
            FW_G2A_uniq_lst=[[FW_G2A_U[u][1],u] for u in Unique_FW_G2A]
            RC_C2T_uniq_lst=[[RC_C2T_U[u][1],u] for u in Unique_RC_C2T]
            RC_G2A_uniq_lst=[[RC_G2A_U[u][1],u] for u in Unique_RC_G2A]
            FW_C2T_uniq_lst.sort()
            RC_C2T_uniq_lst.sort()
            FW_G2A_uniq_lst.sort()
            RC_G2A_uniq_lst.sort()
            FW_C2T_uniq_lst=[x[1] for x in FW_C2T_uniq_lst]
            RC_C2T_uniq_lst=[x[1] for x in RC_C2T_uniq_lst]
            FW_G2A_uniq_lst=[x[1] for x in FW_G2A_uniq_lst]
            RC_G2A_uniq_lst=[x[1] for x in RC_G2A_uniq_lst]

            del Unique_FW_C2T
            del Unique_FW_G2A
            del Unique_RC_C2T
            del Unique_RC_G2A


            #----------------------------------------------------------------
            # Post-filtering reads
            # ---- FW_C2T  ----
            FW_regions = dict()
            gseq = dict()
            chr_length = dict()
            for header in FW_C2T_uniq_lst :
                _, mapped_chr, mapped_location, cigar = FW_C2T_U[header]
                original_BS = original_bs_reads[header]
                if mapped_chr not in FW_regions :
                    FW_regions[mapped_chr] = deserialize(db_d(FWD_MAPPABLE_REGIONS(mapped_chr)))
                if mapped_chr not in gseq :
                    gseq[mapped_chr] = deserialize(db_d(mapped_chr))
                    chr_length[mapped_chr] = len(gseq[mapped_chr])

                r_start, r_end, g_len = get_read_start_end_and_genome_length(cigar)
                all_mapped+=1
                FR = "+FW"
                mapped_strand = "+"
                origin_genome, next2bp, output_genome = get_genomic_sequence(gseq[mapped_chr],
                                                                             mapped_location,
                                                                             mapped_location + g_len,
                                                                             mapped_strand)
                checking_first_C = (output_genome[1:6] == "C_CGG")
                r_aln, g_aln = cigar_to_alignment(cigar, original_BS, origin_genome)

                if len(r_aln) == len(g_aln) :
                    my_region_serial = 0
                    if checking_first_C :
                        my_region_serial, my_region_start, my_region_end = my_mappable_region(FW_regions[mapped_chr],
                                                                                              mapped_location + 1,
                                                                                              FR)
                    if my_region_serial == 0 : # still be 0
                        # for some cases, read has no tags; searching the upstream sequence for CCGG tags
                        # print "[For debug]: FW read has no tags"
                        upstream, _, _ = get_genomic_sequence(gseq[mapped_chr], mapped_location - 500,
                                                              mapped_location + 4, '-')
                        # print "[For debug]: upstream=", upstream
                        # print "[For debug]: mapped_locaiton: %d, CCGG_pos: %d" %(mapped_location, upstream.find('CCGG'))
                        my_region_serial, my_region_start, my_region_end = my_mappable_region(FW_regions[mapped_chr],
                                                                                              mapped_location + 2 - upstream.find('CCGG'), FR)
                        #if my_region_serial == 0 :
                        #    print "[For debug]: chr=", mapped_chr
                        #    print "[For debug]: +FW read still can not find fragment serial"
                        # Tip: sometimes "my_region_serial" is still 0 ...


                    N_mismatch = N_MIS(r_aln, g_aln)
                    if N_mismatch <= int(max_mismatch_no) :
                        all_mapped_passed += 1
                        methy = methy_seq(r_aln, g_aln + next2bp)
                        mC_lst, uC_lst = mcounts(methy, mC_lst, uC_lst)
                        #---XS FILTER----------------
                        XS = 0
                        nCH = methy.count('y') + methy.count('z')
                        nmCH = methy.count('Y') + methy.count('Z')
                        if( (nmCH>XS_count) and nmCH/float(nCH+nmCH)>XS_pct ) :
                            XS = 1
                        num_mapped_FW_C2T += 1
                        outfile.store(header, N_mismatch, FR, mapped_chr, mapped_strand,
                                      mapped_location, cigar, original_BS, methy, XS,
                                      output_genome = output_genome,
                                      rrbs = True,
                                      my_region_serial = my_region_serial,
                                      my_region_start = my_region_start,
                                      my_region_end = my_region_end)
                else :
                    print "[For debug]: reads not in same lengths"


            # ---- RC_C2T ----
            RC_regions = dict()
            for header in RC_C2T_uniq_lst :
                _, mapped_chr, mapped_location, cigar = RC_C2T_U[header]
                original_BS = original_bs_reads[header]
                if mapped_chr not in RC_regions :
                    RC_regions[mapped_chr] = deserialize(db_d(REV_MAPPABLE_REGIONS(mapped_chr)))
                if mapped_chr not in gseq :
                    gseq[mapped_chr] = deserialize(db_d(mapped_chr))
                    chr_length[mapped_chr] = len(gseq[mapped_chr])

                r_start, r_end, g_len = get_read_start_end_and_genome_length(cigar)
                mapped_location = chr_length[mapped_chr] - mapped_location - g_len
                all_mapped+=1
                FR = "-FW"
                mapped_strand = "-"
                origin_genome, next2bp, output_genome = get_genomic_sequence(gseq[mapped_chr],
                                                                             mapped_location,
                                                                             mapped_location + g_len,
                                                                             mapped_strand)
                checking_first_C = (output_genome[1:6] == "C_CGG")
                r_aln, g_aln = cigar_to_alignment(cigar, original_BS, origin_genome)

                if len(r_aln) == len(g_aln) : # and checking_first_C:
                    my_region_serial = 0
                    if checking_first_C :
                        my_region_serial, my_region_start, my_region_end = my_mappable_region(RC_regions[mapped_chr],
                                                                                              mapped_location + g_len + 1,
                                                                                              FR)
                    if my_region_serial == 0 : # still be 0
                        # for some cases, read has no tags; searching the upstream sequence for CCGG tags
                        # print "[For debug]: RC Read has no tags"
                        upstream, _, _ = get_genomic_sequence(gseq[mapped_chr], mapped_location - 3,
                                                              mapped_location + 500, '+')
                        #print "[For debug]: upstream = ", upstream
                        #print "[For debug]: mapped_location = %d, CCGG_pos: %d" % (mapped_location, upstream.find('CCGG'))
                        my_region_serial, my_region_start, my_region_end = my_mappable_region(RC_regions[mapped_chr],
                                                                                              mapped_location + 1 + upstream.find('CCGG'), FR)
                        #if my_region_serial == 0 :
                        #    print "[For debug]: chr=", mapped_chr
                        #    print "[For debug]: -FW read still cannot find fragment serial"


                    N_mismatch = N_MIS(r_aln, g_aln)
                    if N_mismatch <= int(max_mismatch_no) :
                        all_mapped_passed += 1
                        methy = methy_seq(r_aln, g_aln + next2bp)
                        mC_lst, uC_lst = mcounts(methy, mC_lst, uC_lst)
                        #---XS FILTER----------------
                        XS = 0
                        nCH = methy.count('y') + methy.count('z')
                        nmCH = methy.count('Y') + methy.count('Z')
                        if( (nmCH>XS_count) and nmCH/float(nCH+nmCH)>XS_pct ) :
                            XS = 1
                        num_mapped_RC_C2T += 1
                        outfile.store(header, N_mismatch, FR, mapped_chr, mapped_strand,
                                      mapped_location, cigar, original_BS, methy, XS,
                                      output_genome = output_genome,
                                      rrbs = True,
                                      my_region_serial = my_region_serial,
                                      my_region_start = my_region_start,
                                      my_region_end = my_region_end)
                else :
                    print "[For debug]: reads not in same lengths"


            # ---- FW_G2A  ----
            FW_regions = dict()
            gseq = dict()
            chr_length = dict()
            for header in FW_G2A_uniq_lst :
                _, mapped_chr, mapped_location, cigar = FW_G2A_U[header]
                original_BS = original_bs_reads[header]
                if mapped_chr not in FW_regions :
                    FW_regions[mapped_chr] = deserialize(db_d(FWD_MAPPABLE_REGIONS(mapped_chr)))
                if mapped_chr not in gseq :
                    gseq[mapped_chr] = deserialize(db_d(mapped_chr))
                    chr_length[mapped_chr] = len(gseq[mapped_chr])
                cigar = list(reversed(cigar))

                r_start, r_end, g_len = get_read_start_end_and_genome_length(cigar)
                all_mapped+=1
                FR = "-RC"
                mapped_strand = "-"
                origin_genome, next2bp, output_genome = get_genomic_sequence(gseq[mapped_chr],
                                                                             mapped_location,
                                                                             mapped_location + g_len,
                                                                             mapped_strand)
                original_BS = reverse_compl_seq(original_BS)  # for RC reads
                checking_first_C = (output_genome[1:6] == "C_CGG")
                r_aln, g_aln = cigar_to_alignment(cigar, original_BS, origin_genome)

                if len(r_aln) == len(g_aln) :
                    my_region_serial = 0
                    if checking_first_C :
                        my_region_serial, my_region_start, my_region_end = my_mappable_region(FW_regions[mapped_chr],
                                                                                              mapped_location + 1,
                                                                                              FR)
                    if my_region_serial == 0 : # still be 0
                        # for some cases, read has no tags; searching the upstream sequence for CCGG tags
                        # print "[For debug]: FW read has no tags"
                        upstream, _, _ = get_genomic_sequence(gseq[mapped_chr], mapped_location - 500,
                                                              mapped_location + 4, '-')
                        # print "[For debug]: upstream=", upstream
                        # print "[For debug]: mapped_locaiton: %d, CCGG_pos: %d" %(mapped_location, upstream.find('CCGG'))
                        my_region_serial, my_region_start, my_region_end = my_mappable_region(FW_regions[mapped_chr],
                                                                                              mapped_location + 2 - upstream.find('CCGG'), FR)
                        #if my_region_serial == 0 :
                        #    print "[For debug]: chr=", mapped_chr
                        #    print "[For debug]: +FW read still can not find fragment serial"
                        # Tip: sometimes "my_region_serial" is still 0 ...


                    N_mismatch = N_MIS(r_aln, g_aln)
                    if N_mismatch <= int(max_mismatch_no) :
                        all_mapped_passed += 1
                        methy = methy_seq(r_aln, g_aln + next2bp)
                        mC_lst, uC_lst = mcounts(methy, mC_lst, uC_lst)
                        #---XS FILTER----------------
                        XS = 0
                        nCH = methy.count('y') + methy.count('z')
                        nmCH = methy.count('Y') + methy.count('Z')
                        if( (nmCH>XS_count) and nmCH/float(nCH+nmCH)>XS_pct ) :
                            XS = 1
                        num_mapped_FW_G2A += 1
                        outfile.store(header, N_mismatch, FR, mapped_chr, mapped_strand,
                                      mapped_location, cigar, original_BS, methy, XS,
                                      output_genome = output_genome,
                                      rrbs = True,
                                      my_region_serial = my_region_serial,
                                      my_region_start = my_region_start,
                                      my_region_end = my_region_end)
                else :
                    print "[For debug]: reads not in same lengths"


            # ---- RC_G2A ----
            RC_regions = dict()
            for header in RC_G2A_uniq_lst :
                _, mapped_chr, mapped_location, cigar = RC_G2A_U[header]
                original_BS = original_bs_reads[header]
                if mapped_chr not in RC_regions :
                    RC_regions[mapped_chr] = deserialize(db_d(REV_MAPPABLE_REGIONS(mapped_chr)))
                if mapped_chr not in gseq :
                    gseq[mapped_chr] = deserialize(db_d(mapped_chr))
                    chr_length[mapped_chr] = len(gseq[mapped_chr])
                cigar = list(reversed(cigar))

                r_start, r_end, g_len = get_read_start_end_and_genome_length(cigar)
                mapped_location = chr_length[mapped_chr] - mapped_location - g_len
                all_mapped+=1
                FR = "+RC"
                mapped_strand = "+"
                origin_genome, next2bp, output_genome = get_genomic_sequence(gseq[mapped_chr],
                                                                             mapped_location,
                                                                             mapped_location + g_len,
                                                                             mapped_strand)
                original_BS = reverse_compl_seq(original_BS)  # for RC reads
                checking_first_C = (output_genome[1:6] == "C_CGG")
                r_aln, g_aln = cigar_to_alignment(cigar, original_BS, origin_genome)

                if len(r_aln) == len(g_aln) : # and checking_first_C:
                    my_region_serial = 0
                    if checking_first_C :
                        my_region_serial, my_region_start, my_region_end = my_mappable_region(RC_regions[mapped_chr],
                                                                                              mapped_location + g_len + 1,
                                                                                              FR)
                    if my_region_serial == 0 : # still be 0
                        # for some cases, read has no tags; searching the upstream sequence for CCGG tags
                        # print "[For debug]: RC Read has no tags"
                        upstream, _, _ = get_genomic_sequence(gseq[mapped_chr], mapped_location - 3,
                                                              mapped_location + 500, '+')
                        #print "[For debug]: upstream = ", upstream
                        #print "[For debug]: mapped_location = %d, CCGG_pos: %d" % (mapped_location, upstream.find('CCGG'))
                        my_region_serial, my_region_start, my_region_end = my_mappable_region(RC_regions[mapped_chr],
                                                                                              mapped_location + 1 + upstream.find('CCGG'), FR)
                        #if my_region_serial == 0 :
                        #    print "[For debug]: chr=", mapped_chr
                        #    print "[For debug]: -FW read still cannot find fragment serial"


                    N_mismatch = N_MIS(r_aln, g_aln)
                    if N_mismatch <= int(max_mismatch_no) :
                        all_mapped_passed += 1
                        methy = methy_seq(r_aln, g_aln + next2bp)
                        mC_lst, uC_lst = mcounts(methy, mC_lst, uC_lst)
                        #---XS FILTER----------------
                        XS = 0
                        nCH = methy.count('y') + methy.count('z')
                        nmCH = methy.count('Y') + methy.count('Z')
                        if( (nmCH>XS_count) and nmCH/float(nCH+nmCH)>XS_pct ) :
                            XS = 1
                        num_mapped_RC_G2A += 1
                        outfile.store(header, N_mismatch, FR, mapped_chr, mapped_strand,
                                      mapped_location, cigar, original_BS, methy, XS,
                                      output_genome = output_genome,
                                      rrbs = True,
                                      my_region_serial = my_region_serial,
                                      my_region_start = my_region_start,
                                      my_region_end = my_region_end)
                else :
                    print "[For debug]: reads not in same lengths"



            # Finished both FW and RC
            logm("Done: %s (%d) \n" % (read_file, no_my_files))
            print "--> %s (%d) "%(read_file, no_my_files)
            del original_bs_reads
            delete_files(WC2T, CC2T, WG2A, CG2A)



    # End of un-directional library

    delete_files(tmp_path)


    logm("O Number of raw reads: %d "% all_raw_reads)
    if all_raw_reads >0:
        logm("O Number of CGG/TGG tagged reads: %d (%1.3f)"%(all_tagged,float(all_tagged)/all_raw_reads))
        for kk in range(len(n_mytag_lst)):
            logm("O Number of raw reads with %s tag: %d (%1.3f)"%(mytag_lst[kk],n_mytag_lst[mytag_lst[kk]],float(n_mytag_lst[mytag_lst[kk]])/all_raw_reads))
        logm("O Number of CGG/TGG reads having adapter removed: %d "%all_tagged_trimed)
        logm("O Number of unique-hits reads for post-filtering: %d"%all_mapped)

        logm("O ------ %d uniquely aligned reads, passed fragment check, with mismatches <= %s"%(all_mapped_passed, max_mismatch_no))
        logm("O Mappability= %1.4f%%"%(100*float(all_mapped_passed)/all_raw_reads))

        if asktag=="Y": # undiretional
            logm(" ---- %7d FW reads mapped to Watson strand"%(num_mapped_FW_C2T) )
            logm(" ---- %7d RC reads mapped to Watson strand"%(num_mapped_FW_G2A) )
            logm(" ---- %7d FW reads mapped to Crick strand"%(num_mapped_RC_C2T) )
            logm(" ---- %7d RC reads mapped to Crick strand"%(num_mapped_RC_G2A) )
            # the variable name 'num_mapped_RC_G2A' seems not consistent with illustration
            # according to literal meaning
        elif asktag=="N": # directional
            logm(" ---- %7d FW reads mapped to Watson strand"%(num_mapped_FW_C2T) )
            logm(" ---- %7d FW reads mapped to Crick strand"%(num_mapped_RC_C2T) )

        n_CG=mC_lst[0]+uC_lst[0]
        n_CHG=mC_lst[1]+uC_lst[1]
        n_CHH=mC_lst[2]+uC_lst[2]

        logm("----------------------------------------------")
        logm("M Methylated C in mapped reads ")
        logm("M mCG %1.3f%%"%((100*float(mC_lst[0])/n_CG) if n_CG != 0 else 0))
        logm("M mCHG %1.3f%%"%((100*float(mC_lst[1])/n_CHG) if n_CHG != 0 else 0))
        logm("M mCHH %1.3f%%"%((100*float(mC_lst[2])/n_CHH) if n_CHH != 0 else 0))
    logm("----------------------------------------------")
    logm("------------------- END ----------------------")
    elapsed(main_read_file)

    close_log()


