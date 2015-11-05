import fileinput, os, time, random, math
from bs_utils.utils import *
from bs_align_utils import *
import gzip

#----------------------------------------------------------------
# Read from the mapped results, return lists of unique / multiple-hit reads
# The function suppose at most 2 hits will be reported in single file
def extract_mapping(ali_file):
    unique_hits = {}
    non_unique_hits = {}

    header0 = ""
    lst = []

    for header, chr, location, no_mismatch, cigar in process_aligner_output(ali_file):
        #------------------------------
        if header != header0:
            #---------- output -----------
            if len(lst) == 1:
                unique_hits[header0] = lst[0]      # [no_mismatch, chr, location]
            elif len(lst) > 1:
                min_lst = min(lst, key = lambda x: x[0])
                max_lst = max(lst, key = lambda x: x[0])

                if min_lst[0] < max_lst[0]:
                    unique_hits[header0] = min_lst
                else:
                    non_unique_hits[header0] = min_lst[0]
                    #print "multiple hit", header, chr, location, no_mismatch, cigar # test
            header0 = header
            lst = [(no_mismatch, chr, location, cigar)]
        else: # header == header0, same header (read id)
            lst.append((no_mismatch, chr, location, cigar))

    if len(lst) == 1:
        unique_hits[header0] = lst[0]      # [no_mismatch, chr, location]
    elif len(lst) > 1:
        min_lst = min(lst, key = lambda x: x[0])
        max_lst = max(lst, key = lambda x: x[0])

        if min_lst[0] < max_lst[0]:
            unique_hits[header0] = min_lst
        else:
            non_unique_hits[header0] = min_lst[0]


    return unique_hits, non_unique_hits


def bs_single_end(main_read_file, asktag, adapter_file, cut1, cut2, no_small_lines,
                  max_mismatch_no, aligner_command, db_path, tmp_path, outfile,
                  XS_pct, XS_count, adapter_mismatch, show_multiple_hit, show_unmapped_hit):
    logm("----------------------------------------------" )
    logm("Read filename: %s" % main_read_file)
    logm("The first base (for mapping): %d" % cut1  )
    logm("The  last base (for mapping): %d" % cut2  )
    logm("Path for short reads aligner: %s" % aligner_command + '\n')
    logm("Reference genome library path: %s" % db_path  )
    if asktag == "Y" :
        logm("Un-directional library" )
    else :
        logm("Directional library")
    # end-of-if
    logm("Number of mismatches allowed: %s" % str(max_mismatch_no)  )
    # adapter : strand-specific or not
    adapter = ""
    adapter_fw = ""
    adapter_rc = ""
    if adapter_file != "":
        try :
            adapter_inf = open(adapter_file, "r")
            if asktag == "N": #<--- directional library
                adapter = adapter_inf.readline()
                adapter_inf.close()
                adapter = adapter.rstrip("\n")[0:10]
            elif asktag == "Y":#<--- un-directional library
                adapter_fw = adapter_inf.readline()
                adapter_rc = adapter_inf.readline()
                adapter_inf.close()
                adapter_fw = adapter_fw.rstrip("\n")[0:10]
                adapter_rc = adapter_rc.rstrip("\n")[-10::]
                if adapter_rc == "" :
                    adapter_rc = reverse_compl_seq(adapter_fw)
            adapter_inf.close()
        except IOError:
            print "[Error] Cannot open adapter file : %s" % adapter_file
            exit(-1)
    if adapter_file != "":
        if asktag == "N": #<--- directional library
            logm("Adapter sequence: %s" % adapter)
        elif asktag == "Y":
            logm("3\' end adapter sequence: %s" % adapter_fw)
            logm("5\' end adapter sequence: %s" % adapter_rc)
    logm("-------------------------------- " )

    # helper method to join fname with tmp_path
    tmp_d = lambda fname: os.path.join(tmp_path, fname)
    db_d = lambda fname:  os.path.join(db_path, fname)

    # splitting the big read file
    input_fname = os.path.split(main_read_file)[1]

    #---- Stats ------------------------------------------------------------
    all_raw_reads = 0
    all_trimmed = 0
    all_mapped = 0
    all_mapped_passed = 0
    all_base_before_trim = 0
    all_base_after_trim = 0
    all_base_mapped = 0
    numbers_premapped_lst = [0, 0, 0, 0]
    numbers_mapped_lst = [0, 0, 0, 0]
    mC_lst = [0, 0, 0]
    uC_lst = [0, 0, 0]
    no_my_files = 0
    #----------------------------------------------------------------

    if show_multiple_hit is not None:
        outf_MH=open(show_multiple_hit,'w')
    if show_unmapped_hit is not None :
        outf_UH=open(show_unmapped_hit,'w')

    for read_file in isplit_file(main_read_file, tmp_d(input_fname)+'-s-', no_small_lines):
#    for read_file in my_files:
        original_bs_reads = {}
        no_my_files+=1
        random_id = ".tmp-"+str(random.randint(1000000, 9999999))
        #-------------------------------------------------------------------
        # un-directional sequencing
        #-------------------------------------------------------------------
        if asktag=="Y":  

            #----------------------------------------------------------------
            outfile2=tmp_d('Trimmed_C2T.fa'+random_id)
            outfile3=tmp_d('Trimmed_G2A.fa'+random_id)

            outf2=open(outfile2,'w')
            outf3=open(outfile3,'w')

            #----------------------------------------------------------------
            # detect format of input file
            try :
                if read_file.endswith(".gz") : # support input file ending with ".gz"
                    read_inf = gzip.open(read_file, "rb")
                else :
                    read_inf=open(read_file,"r")
            except IOError :
                print "[Error] Cannot open input file : %s" % read_file
                exit(-1)

            logm("Start reading and trimming the input sequences")
            oneline = read_inf.readline()
            if oneline == "" :
                oneline = "NNNN"
            l = oneline.split()
            input_format = ""
            if oneline[0]=="@":
                input_format = "fastq"
            elif len(l)==1 and oneline[0]!=">":
                input_format = "seq"
            elif len(l)==11:
                input_format = "qseq"
            elif oneline[0]==">":
                input_format = "fasta"
            read_inf.close()

            #----------------------------------------------------------------
            # read sequence, remove adapter and convert
            read_id = ""
            seq = ""
            seq_ready = "N"
            line_no = 0
            fw_trimmed = 0
            rc_trimmed = 0
            for line in fileinput.input(read_file, openhook=fileinput.hook_compressed): # allow input with .gz
                if line == "" : # fix bug for empty input line
                    line = "NNNN"
                l = line.split()
                line_no += 1
                if input_format=="seq":
                    all_raw_reads += 1
                    read_id = str(all_raw_reads)
                    read_id = read_id.zfill(12)
                    seq = l[0]
                    seq_ready = "Y"
                elif input_format=="fastq":
                    l_fastq = math.fmod(line_no, 4)
                    if l_fastq == 1 :
                        all_raw_reads += 1
                        read_id = l[0][1:]
                        seq_ready = "N"
                    elif l_fastq == 2 :
                        seq = l[0]
                        seq_ready = "Y"
                    else :
                        seq = ""
                        seq_ready = "N"
                elif input_format=="qseq":
                    all_raw_reads += 1
                    read_id = str(all_raw_reads)
                    read_id = read_id.zfill(12)
                    seq = l[8]
                    seq_ready = "Y"
                elif input_format=="fasta" :
                    l_fasta = math.fmod(line_no,2)
                    if l_fasta==1:
                        all_raw_reads += 1
                        read_id = l[0][1:]
                        seq = ""
                        seq_ready = "N"
                    elif l_fasta==0 :
                        seq = l[0]
                        seq_ready = "Y"

                #----------------------------------------------------------------
                if seq_ready=="Y":
                    seq=seq[cut1-1:cut2] #<---- selecting 0..52 from 1..72  -e 52
                    seq=seq.upper()
                    seq=seq.replace(".","N")

                    # striping BS adapter from 3' read
                    all_base_before_trim += len(seq)
                    if (adapter_fw !="") or (adapter_rc !="") :
                        new_read = RemoveAdapter(seq, adapter_fw, adapter_mismatch)
                        if len(new_read) < len(seq) :
                            fw_trimmed += 1
                        new_read_len = len(new_read)
                        #print new_read
                        new_read = Remove_5end_Adapter(new_read, adapter_rc, adapter_mismatch)
                        new_read = RemoveAdapter(new_read, adapter_fw, adapter_mismatch)
                        if len(new_read) < new_read_len :
                            rc_trimmed += 1
                        #print new_read
                        if len(new_read) < len(seq) :
                            all_trimmed += 1
                        seq = new_read
                    all_base_after_trim += len(seq)
                    if len(seq)<=4:
                        seq=''.join(["N" for x in xrange(cut2-cut1+1)])

                    #---------  trimmed_raw_BS_read  ------------------
                    original_bs_reads[read_id] = seq

                    #---------  FW_C2T  ------------------
                    outf2.write('>%s\n%s\n' % (read_id, seq.replace("C","T")))
                    #---------  RC_G2A  ------------------
                    outf3.write('>%s\n%s\n' % (read_id, seq.replace("G","A")))

            fileinput.close()

            outf2.close()
            outf3.close()

            delete_files(read_file)
            logm("Reads trimmed from 3\' end : %d " % fw_trimmed)
            logm("Reads trimmed from 5\' end : %d " % rc_trimmed)
           #--------------------------------------------------------------------------------
            # Bowtie mapping
            #-------------------------------------------------------------------------------
            logm("Start mapping")
            WC2T=tmp_d("W_C2T_m"+str(max_mismatch_no)+".mapping"+random_id)
            CC2T=tmp_d("C_C2T_m"+str(max_mismatch_no)+".mapping"+random_id)
            WG2A=tmp_d("W_G2A_m"+str(max_mismatch_no)+".mapping"+random_id)
            CG2A=tmp_d("C_G2A_m"+str(max_mismatch_no)+".mapping"+random_id)

        #    print aligner_command % {'int_no_mismatches' : int_no_mismatches,
        #                             'reference_genome' : os.path.join(db_path,'W_C2T'),
        #                             'input_file' : outfile2,
        #                             'output_file' : WC2T}

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


            delete_files(outfile2, outfile3)


            #--------------------------------------------------------------------------------
            # Post processing
            #--------------------------------------------------------------------------------

            FW_C2T_U,FW_C2T_R=extract_mapping(WC2T)
            RC_G2A_U,RC_G2A_R=extract_mapping(CG2A)

            FW_G2A_U,FW_G2A_R=extract_mapping(WG2A)
            RC_C2T_U,RC_C2T_R=extract_mapping(CC2T)

            #----------------------------------------------------------------
            # get unique-hit reads
            #----------------------------------------------------------------
            Union_set=set(FW_C2T_U.iterkeys()) | set(RC_G2A_U.iterkeys()) | set(FW_G2A_U.iterkeys()) | set(RC_C2T_U.iterkeys())

            Unique_FW_C2T=set() # +
            Unique_RC_G2A=set() # +
            Unique_FW_G2A=set() # -
            Unique_RC_C2T=set() # -
            Multiple_hits=set()


            for x in Union_set:
                _list=[]
                for d in [FW_C2T_U, RC_G2A_U, FW_G2A_U, RC_C2T_U]:
                    mis_lst=d.get(x,[99])
                    mis=int(mis_lst[0])
                    _list.append(mis)
                for d in [FW_C2T_R, RC_G2A_R, FW_G2A_R, RC_C2T_R]:
                    mis=d.get(x,99) 
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
                    # if mini_index = 4,5,6,7, indicating multiple hits
                    else :
                        Multiple_hits.add(x)
                else :
                    Multiple_hits.add(x)
            # write reads rejected by Multiple Hits to file
            if show_multiple_hit is not None :
                #outf_MH=open(show_multiple_hit,'w')
                for i in Multiple_hits :
                    outf_MH.write(">%s\n" % i)
                    outf_MH.write("%s\n" % original_bs_reads[i])
                #outf_MH.close()

            # write unmapped reads to file
            if show_unmapped_hit is not None :
                #outf_UH=open(show_unmapped_hit,'w')
                for i in original_bs_reads :
                    if i not in Union_set :
                        outf_UH.write(">%s\n" % i)
                        outf_UH.write("%s\n" % original_bs_reads[i])
                #outf_UH.close()

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
            #----------------------------------------------------------------

            numbers_premapped_lst[0] += len(Unique_FW_C2T)
            numbers_premapped_lst[1] += len(Unique_RC_G2A)
            numbers_premapped_lst[2] += len(Unique_FW_G2A)
            numbers_premapped_lst[3] += len(Unique_RC_C2T)

            del Unique_FW_C2T
            del Unique_FW_G2A
            del Unique_RC_C2T
            del Unique_RC_G2A


            #----------------------------------------------------------------

            nn=0
            gseq = dict()
            chr_length = dict()
            for ali_unique_lst, ali_dic in [(FW_C2T_uniq_lst,FW_C2T_U),
                                            (RC_G2A_uniq_lst,RC_G2A_U),
                                            (FW_G2A_uniq_lst,FW_G2A_U),
                                            (RC_C2T_uniq_lst,RC_C2T_U)]:
                nn += 1

                for header in ali_unique_lst:

                    _, mapped_chr, mapped_location, cigar = ali_dic[header]

                    original_BS = original_bs_reads[header]
                    #-------------------------------------
                    if mapped_chr not in gseq:
                        gseq[mapped_chr] =  deserialize(db_d(mapped_chr))
                        chr_length[mapped_chr] = len(gseq[mapped_chr])

                    if nn == 2 or nn == 3:
                        cigar = list(reversed(cigar))
                    r_start, r_end, g_len = get_read_start_end_and_genome_length(cigar)


                    all_mapped += 1

                    if nn == 1: # +FW mapped to + strand:
                        FR = "+FW"
                        mapped_strand="+"

                    elif nn == 2:  # +RC mapped to + strand:
                        FR = "+RC" # RC reads from -RC reflecting the methylation status on Watson strand (+)
                        mapped_location = chr_length[mapped_chr] - mapped_location - g_len
                        mapped_strand = "+"
                        original_BS = reverse_compl_seq(original_BS)  # for RC reads

                    elif nn == 3:  						# -RC mapped to - strand:
                        mapped_strand = "-"
                        FR = "-RC" # RC reads from +RC reflecting the methylation status on Crick strand (-)
                        original_BS = reverse_compl_seq(original_BS)  # for RC reads

                    elif nn == 4: 						# -FW mapped to - strand:
                        mapped_strand = "-"
                        FR = "-FW"
                        mapped_location = chr_length[mapped_chr] - mapped_location - g_len

                    origin_genome, next, output_genome = get_genomic_sequence(gseq[mapped_chr], mapped_location, mapped_location + g_len, mapped_strand)

                    r_aln, g_aln = cigar_to_alignment(cigar, original_BS, origin_genome)


                    if len(r_aln)==len(g_aln):
                        N_mismatch = N_MIS(r_aln, g_aln)
#                        if N_mismatch <= int(max_mismatch_no):
                        mm_no=float(max_mismatch_no)
                        if (mm_no>=1 and N_mismatch<=mm_no) or (mm_no<1 and N_mismatch<=(mm_no*len(r_aln)) ):
                            numbers_mapped_lst[nn-1] += 1
                            all_mapped_passed += 1
                            methy = methy_seq(r_aln, g_aln + next)
                            mC_lst, uC_lst = mcounts(methy, mC_lst, uC_lst)

                            #---XS FILTER----------------
                            XS = 0
                            nCH = methy.count('y') + methy.count('z')
                            nmCH = methy.count('Y') + methy.count('Z')
                            if( (nmCH>XS_count) and nmCH/float(nCH+nmCH)>XS_pct ) :
                                XS = 1

                            outfile.store(header, N_mismatch, FR, mapped_chr, mapped_strand, mapped_location, cigar, original_BS, methy, XS, output_genome = output_genome)
                            all_base_mapped += len(original_BS)

            #----------------------------------------------------------------
            logm("--> %s (%d) "%(read_file, no_my_files))
            delete_files(WC2T, WG2A, CC2T, CG2A)



        #--------------------------------------------------------------------
        # directional sequencing
        #--------------------------------------------------------------------

        if asktag=="N":  
            #----------------------------------------------------------------
            outfile2=tmp_d('Trimmed_C2T.fa'+random_id)
            outf2=open(outfile2,'w')
            #----------------------------------------------------------------
            try :
                if read_file.endswith(".gz") : # support input file ending with ".gz"
                    read_inf = gzip.open(read_file, "rb")
                else :
                    read_inf=open(read_file,"r")
            except IOError :
                print "[Error] Cannot open input file : %s" % read_file
                exit(-1)

            logm("Start reading and trimming the input sequences")
            oneline = read_inf.readline()
            if oneline == "" :
                oneline = "NNNN"
            l = oneline.split()
            input_format = ""
            if oneline[0]=="@":
                input_format = "fastq"
            elif len(l)==1 and oneline[0]!=">":
                input_format = "seq"
            elif len(l)==11:
                input_format = "qseq"
            elif oneline[0]==">":
                input_format = "fasta"
            read_inf.close()

            #print "detected data format: %s"%(input_format)
            #----------------------------------------------------------------
            read_id=""
            seq=""
            seq_ready="N"
            line_no = 0
            for line in fileinput.input(read_file, openhook=fileinput.hook_compressed):
                if l == "" :
                    l = "NNNN"
                l = line.split()
                line_no += 1
                if input_format=="seq":
                    all_raw_reads += 1
                    read_id = str(all_raw_reads)
                    read_id = read_id.zfill(12)
                    seq = l[0]
                    seq_ready = "Y"
                elif input_format=="fastq":
                    l_fastq = math.fmod(line_no, 4)
                    if l_fastq == 1 :
                        all_raw_reads += 1
                        read_id = l[0][1:]
                        seq_ready = "N"
                    elif l_fastq == 2 :
                        seq = l[0]
                        seq_ready = "Y"
                    else :
                        seq = ""
                        seq_ready = "N"
                elif input_format=="qseq":
                    all_raw_reads += 1
                    read_id = str(all_raw_reads)
                    read_id = read_id.zfill(12)
                    seq = l[8]
                    seq_ready = "Y"
                elif input_format=="fasta" :
                    l_fasta = math.fmod(line_no,2)
                    if l_fasta==1:
                        all_raw_reads += 1
                        read_id = l[0][1:]
                        seq = ""
                        seq_ready = "N"
                    elif l_fasta==0 :
                        seq = l[0]
                        seq_ready = "Y"
                #--------------------------------
                if seq_ready=="Y":
                    seq=seq[cut1-1:cut2] #<---selecting 0..52 from 1..72  -e 52
                    seq=seq.upper()
                    seq=seq.replace(".","N")

                    #--striping adapter from 3' read -------
                    all_base_before_trim += len(seq)
                    if adapter != "":
                        new_read = RemoveAdapter(seq, adapter, adapter_mismatch)
                        if len(new_read) < len(seq) :
                            all_trimmed += 1
                        seq = new_read
                    all_base_after_trim += len(seq)
                    if len(seq)<=4:
                        seq = "N" * (cut2-cut1+1)

                    #---------  trimmed_raw_BS_read  ------------------
                    original_bs_reads[read_id] = seq


                    #---------  FW_C2T  ------------------
                    outf2.write('>%s\n%s\n' % (read_id, seq.replace("C","T")))

            fileinput.close()

            outf2.close()
            delete_files(read_file)

            #--------------------------------------------------------------------------------
            # Bowtie mapping
            #--------------------------------------------------------------------------------
            logm("Start mapping")
            WC2T=tmp_d("W_C2T_m"+str(max_mismatch_no)+".mapping"+random_id)
            CC2T=tmp_d("C_C2T_m"+str(max_mismatch_no)+".mapping"+random_id)

            run_in_parallel([ aligner_command % {'reference_genome' : os.path.join(db_path,'W_C2T'),
                                                  'input_file' : outfile2,
                                                  'output_file' : WC2T},
                              aligner_command % {'reference_genome' : os.path.join(db_path,'C_C2T'),
                                                  'input_file' : outfile2,
                                                  'output_file' : CC2T} ])

            delete_files(outfile2)

            #--------------------------------------------------------------------------------
            # Post processing
            #--------------------------------------------------------------------------------


            FW_C2T_U, FW_C2T_R = extract_mapping(WC2T)
            RC_C2T_U, RC_C2T_R = extract_mapping(CC2T)

            #----------------------------------------------------------------
            # get uniq-hit reads
            #----------------------------------------------------------------
            Union_set = set(FW_C2T_U.iterkeys()) | set(RC_C2T_U.iterkeys())

            Unique_FW_C2T = set() # +
            Unique_RC_C2T = set() # -
            Multiple_hits=set()
            # write reads rejected by Multiple Hits to file

            for x in Union_set:
                _list=[]
                for d in [FW_C2T_U,RC_C2T_U]:
                    mis_lst=d.get(x,[99])
                    mis=int(mis_lst[0])
                    _list.append(mis)
                for d in [FW_C2T_R,RC_C2T_R]:
                    mis=d.get(x,99)
                    _list.append(mis)
                mini=min(_list)
                #print _list
                if _list.count(mini)==1:
                    mini_index=_list.index(mini)
                    if mini_index==0:
                        Unique_FW_C2T.add(x)
                    elif mini_index==1:
                        Unique_RC_C2T.add(x)
                    else:
                        Multiple_hits.add(x)
                else :
                    Multiple_hits.add(x)
            # write reads rejected by Multiple Hits to file
            if show_multiple_hit is not None:
                #outf_MH=open(show_multiple_hit,'w')
                for i in Multiple_hits :
                    outf_MH.write(">%s\n" % i)
                    outf_MH.write("%s\n" % original_bs_reads[i])
                #outf_MH.close()

            # write unmapped reads to file
            if show_unmapped_hit is not None :
                #outf_UH=open(show_unmapped_hit,'w')
                for i in original_bs_reads :
                    if i not in Union_set :
                        outf_UH.write(">%s\n" % i)
                        outf_UH.write("%s\n" % original_bs_reads[i])
                #outf_UH.close()

            FW_C2T_uniq_lst=[[FW_C2T_U[u][1],u] for u in Unique_FW_C2T]
            RC_C2T_uniq_lst=[[RC_C2T_U[u][1],u] for u in Unique_RC_C2T]
            FW_C2T_uniq_lst.sort()
            RC_C2T_uniq_lst.sort()
            FW_C2T_uniq_lst=[x[1] for x in FW_C2T_uniq_lst]
            RC_C2T_uniq_lst=[x[1] for x in RC_C2T_uniq_lst]


            #----------------------------------------------------------------

            numbers_premapped_lst[0] += len(Unique_FW_C2T)
            numbers_premapped_lst[1] += len(Unique_RC_C2T)

            #----------------------------------------------------------------

            nn = 0
            gseq = dict()
            chr_length = dict()
            for ali_unique_lst, ali_dic in [(FW_C2T_uniq_lst,FW_C2T_U),(RC_C2T_uniq_lst,RC_C2T_U)]:
                nn += 1
                for header in ali_unique_lst:
                    _, mapped_chr, mapped_location, cigar = ali_dic[header]
                    original_BS = original_bs_reads[header]
                    #-------------------------------------
                    if mapped_chr not in gseq :
                        gseq[mapped_chr] = deserialize(db_d(mapped_chr))
                        chr_length[mapped_chr] = len(gseq[mapped_chr])

                    r_start, r_end, g_len = get_read_start_end_and_genome_length(cigar)

                    all_mapped+=1
                    if nn == 1: 	# +FW mapped to + strand:
                        FR = "+FW"
                        mapped_strand = "+"
                    elif nn == 2: 	# -FW mapped to - strand:
                        mapped_strand = "-"
                        FR = "-FW"
                        mapped_location = chr_length[mapped_chr] - mapped_location - g_len


                    origin_genome, next, output_genome = get_genomic_sequence(gseq[mapped_chr], mapped_location, mapped_location + g_len, mapped_strand)
                    r_aln, g_aln = cigar_to_alignment(cigar, original_BS, origin_genome)

                    if len(r_aln) == len(g_aln):
                        N_mismatch = N_MIS(r_aln, g_aln) #+ original_BS_length - (r_end - r_start) # mismatches in the alignment + soft clipped nucleotides
                        mm_no=float(max_mismatch_no)
                        if (mm_no>=1 and N_mismatch<=mm_no) or (mm_no<1 and N_mismatch<=(mm_no*len(r_aln)) ):
                            numbers_mapped_lst[nn-1] += 1
                            all_mapped_passed += 1
                            methy = methy_seq(r_aln, g_aln+next)
                            mC_lst, uC_lst = mcounts(methy, mC_lst, uC_lst)

                            #---XS FILTER----------------
                            XS = 0
                            nCH = methy.count('y') + methy.count('z')
                            nmCH = methy.count('Y') + methy.count('Z')
                            if( (nmCH>XS_count) and nmCH/float(nCH+nmCH)>XS_pct ) :
                                XS = 1

                            outfile.store(header, N_mismatch, FR, mapped_chr, mapped_strand, mapped_location, cigar, original_BS, methy, XS, output_genome = output_genome)
                            all_base_mapped += len(original_BS)

            #----------------------------------------------------------------
            logm("--> %s (%d) "%(read_file,no_my_files))
            delete_files(WC2T, CC2T)


    #----------------------------------------------------------------

    delete_files(tmp_path)

    if show_multiple_hit is not None:
        outf_MH.close()

    if show_unmapped_hit is not None :
        outf_UH.close()

    logm("----------------------------------------------" )
    logm("Number of raw reads: %d" % all_raw_reads)
    if all_raw_reads > 0 :
        logm("Number of bases in total: %d " % all_base_before_trim)
        if (asktag == "N" and adapter != "") or (asktag == "Y" and adapter_fw != "") :
            logm("Number of reads having adapter removed: %d" % all_trimmed )
            logm("Number of bases after trimming the adapters: %d (%1.3f)" % (all_base_after_trim, float(all_base_after_trim)/all_base_before_trim) )
        logm("Number of reads are rejected because of multiple hits: %d" % len(Multiple_hits) )
        logm("Number of unique-hits reads (before post-filtering): %d" % all_mapped)
        if asktag == "Y":
            logm("  %7d FW reads mapped to Watson strand (before post-filtering)" % (numbers_premapped_lst[0]) )
            logm("  %7d RC reads mapped to Watson strand (before post-filtering)" % (numbers_premapped_lst[1]) )
            logm("  %7d FW reads mapped to Crick strand (before post-filtering)" % (numbers_premapped_lst[2]) )
            logm("  %7d RC reads mapped to Crick strand (before post-filtering)" % (numbers_premapped_lst[3]) )
        elif asktag == "N":
            logm("  %7d FW reads mapped to Watson strand (before post-filtering)" % (numbers_premapped_lst[0]) )
            logm("  %7d FW reads mapped to Crick strand (before post-filtering)" % (numbers_premapped_lst[1]) )

        logm("Post-filtering %d uniquely aligned reads with mismatches <= %s" % (all_mapped_passed, max_mismatch_no) )
        if asktag == "Y":
            logm("  %7d FW reads mapped to Watson strand" % (numbers_mapped_lst[0]) )
            logm("  %7d RC reads mapped to Watson strand" % (numbers_mapped_lst[1]) )
            logm("  %7d FW reads mapped to Crick strand" % (numbers_mapped_lst[2]) )
            logm("  %7d RC reads mapped to Crick strand" % (numbers_mapped_lst[3]) )
        elif asktag == "N":
            logm("  %7d FW reads mapped to Watson strand" % (numbers_mapped_lst[0]) )
            logm("  %7d FW reads mapped to Crick strand" % (numbers_mapped_lst[1]) )
        logm("Mappability = %1.4f%%" % (100*float(all_mapped_passed)/all_raw_reads) )
        logm("Total bases of uniquely mapped reads : %7d" % all_base_mapped )
        #
        n_CG  = mC_lst[0] + uC_lst[0]
        n_CHG = mC_lst[1] + uC_lst[1]
        n_CHH = mC_lst[2] + uC_lst[2]
        #
        logm("----------------------------------------------" )
        logm("Methylated C in mapped reads ")
        #
        logm(" mCG  %1.3f%%" % ((100*float(mC_lst[0])/n_CG) if n_CG != 0 else 0))
        logm(" mCHG %1.3f%%" % ((100*float(mC_lst[1])/n_CHG) if n_CHG != 0 else 0))
        logm(" mCHH %1.3f%%" % ((100*float(mC_lst[2])/n_CHH) if n_CHH != 0 else 0))
        #
    logm("----------------------------------------------" )
    logm("File : %s" % main_read_file )
    elapsed("Resource / CPU time")
    logm("------------------- END --------------------" )
    close_log()



