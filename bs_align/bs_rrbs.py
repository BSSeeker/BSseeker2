import fileinput, random, math, os.path
from bs_index.rrbs_build import FWD_MAPPABLE_REGIONS, REV_MAPPABLE_REGIONS
from bs_utils.utils import *

from bs_align.bs_single_end import extract_mapping
from bs_align_utils import *

def my_mappable_region(chr_regions, mapped_location, FR): # start_position (first C), end_position (last G), serial, sequence
    #print len(chr_regions)
    out_serial = 0
    out_start = -1
    out_end = -1
    #print "mapped_location:", mapped_location
    if FR == "+FW" or FR == "-RC":
        my_location = str(mapped_location)
        if my_location in chr_regions:
            my_lst = chr_regions[my_location]
            out_start = int(my_location)
            out_end = my_lst[0]
            out_serial = my_lst[1]
        #else :
        #    print "[For debug]: +FW location %s cannot be found" % my_location
    elif FR == "-FW" or FR == "+RC":
        my_location = str(mapped_location)
        if my_location in chr_regions:
            my_lst = chr_regions[my_location]
            out_end = int(my_location)
            out_start = my_lst[0]
            out_serial = my_lst[1]
        #else :
        #    print "[For debug]: -FW location %s cannot be found" % my_location

    return out_serial, out_start, out_end


#----------------------------------------------------------------

def bs_rrbs(main_read_file, asktag, adapter_file, cut_s, cut_e, no_small_lines, max_mismatch_no,
            aligner_command, db_path, tmp_path, outfile, XS_pct, XS_count, XSteve, adapter_mismatch,
            show_multiple_hit, show_unmapped_hit, cut_format="C-CGG"
            ):
    #----------------------------------------------------------------
    # For double enzyme: cut_format="C-CGG,A-CTG"; ApekI:"G^CWGC"
    #cut_context = re.sub("-", "", cut_format)
    # Ex. cut_format="C-CGG,AT-CG,G-CWGC"
    cut_format_lst = EnumerateIUPAC(cut_format.upper().split(",")) # ['G-CAGC', 'AT-CG', 'C-CGG', 'G-CTGC']
    cut_context = [i.replace("-","") for i in cut_format_lst] # ['GCAGC', 'ATCG', 'CCGG', 'GCTGC']
    cut5_context = [re.match( r'(.*)\-(.*)', i).group(1) for i in cut_format_lst] # ['G', 'AT', 'C', 'G']
    cut3_context = [re.match( r'(.*)\-(.*)', i).group(2) for i in cut_format_lst] # ['CAGC', 'CG', 'CGG', 'CTGC']
    cut_len = [len(i) for i in cut_context] # [5, 4, 4, 5]
    min_cut5_len = min([len(i) for i in cut5_context])
    #print cut_format_lst
    #print cut_format
    #print cut5_context
    cut_tag_lst = Enumerate_C_to_CT(cut_format_lst) # ['G-TTGC', 'AT-TG', 'G-CAGT', 'T-CGG', 'G-TAGC', 'C-TGG', 'G-CAGC', 'G-CTGC', 'AT-CG', 'T-TGG', 'G-TTGT', 'G-TAGT', 'C-CGG', 'G-CTGT']
    cut5_tag_lst = [re.match(r'(.*)\-(.*)', i).group(1) for i in cut_tag_lst]
    cut3_tag_lst = [re.match(r'(.*)\-(.*)', i).group(2) for i in cut_tag_lst]
    check_pattern = [ i[-2:]+"_"+j for i,j in zip(cut5_tag_lst, cut3_tag_lst) ]
    #print "======="
    #print cut_tag_lst
    #print cut3_tag_lst
    #print cut5_tag_lst
    #print check_pattern
    # set region[gx,gy] for checking_genome_context
    gx = [ 0 if j>2 else 2-j for j in [len(i) for i in cut5_tag_lst] ] # [XC-CGG]
    gy = [ 3+len(i) for i in cut3_tag_lst ]
    #----------------------------------------------------------------
    # helper method to join fname with tmp_path
    tmp_d = lambda fname: os.path.join(tmp_path, fname)
    db_d = lambda fname: os.path.join(db_path, fname)
    MAX_TRY = 500 # For finding the serial_no
    logm("----------------------------------------------" )
    adapter = ""
    adapter_fw = ""
    adapter_rc = ""
    if adapter_file:
        try :
            adapter_inf = open(adapter_file,"r")
            if asktag == "N": #<--- directional library
                adapter = adapter_inf.readline().strip()
                adapter = adapter[0:10] # only use first 10bp of adapter
                adapter_inf.close()
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
            print "[Error] Cannot find adapter file : %s !" % adapter_file
            exit(-1)
    # end-of-if
    logm("Read filename: %s" % main_read_file)
    logm("The first base (for mapping): %d" % cut_s  )
    logm("The  last base (for mapping): %d" % cut_e  )
    logm("Path for short reads aligner: %s" % aligner_command + '\n')
    logm("Reference genome library path: %s" % db_path  )
    if asktag == "Y" :
        logm("Un-directional library" )
    else :
        logm("Directional library")
    # end-of-if
    logm("Number of mismatches allowed: %s" % max_mismatch_no  )

    if adapter_file != "":
        if asktag == "N": #<--- directional library
            logm("Adapter sequence: %s" % adapter)
        elif asktag == "Y":
            logm("3\' end adapter sequence: %s" % adapter_fw)
            logm("5\' end adapter sequence: %s" % adapter_rc)
    #if adapter_file !="":
    #    logm("Adapter seq: %s" % adapter)
    logm("-------------------------------- " )

    all_raw_reads=0
    all_tagged=0
    all_tagged_trimmed=0
    all_base_before_trim=0
    all_base_after_trim=0
    all_base_mapped=0
    all_mapped=0
    all_mapped_passed=0
    n_cut_tag_lst={}
    #print cut3_tag_lst
    for x in cut3_tag_lst:
        n_cut_tag_lst[x]=0

    mC_lst=[0,0,0]
    uC_lst=[0,0,0]

    no_my_files=0

    num_mapped_FW_C2T = 0
    num_mapped_RC_C2T = 0
    num_mapped_FW_G2A = 0
    num_mapped_RC_G2A = 0

    # Count of nucleotides, which should be cut before the adapters
    Extra_base_cut_5end_adapter = max([ abs(len(i)-len(j)) for i,j in zip(cut5_context, cut3_context)])

    #===============================================
    # directional sequencing
    #===============================================

    if show_multiple_hit is not None :
        outf_MH=open(show_multiple_hit,'w')

    if show_unmapped_hit is not None :
        outf_UH=open(show_unmapped_hit,'w')

    if asktag=="N" :
        #----------------------------------------------------------------
        logm("Start reading and trimming the input sequences")

        input_fname = os.path.split(main_read_file)[1]
        for read_file in isplit_file(main_read_file, tmp_d(input_fname)+'-s-', no_small_lines):

            logm("Processing read file: %s" % read_file)
            original_bs_reads = {}
            no_my_files+=1
            random_id = ".tmp-"+str(random.randint(1000000,9999999))
            outfile2=tmp_d('Trimmed_C2T.fa'+random_id)

            outf2=open(outfile2,'w')

            #--- Checking input format ------------------------------------------
            try :
                if read_file.endswith(".gz") : # support input file ending with ".gz"
                    read_inf = gzip.open(read_file, "rb")
                else :
                    read_inf=open(read_file,"r")
            except IOError:
                print "[Error] Cannot open input file : %s" % read_file
                exit(-1)

            oneline = read_inf.readline()
            if oneline == "" :
                oneline = "NNNN"
            l = oneline.split()
            input_format = ""
            if oneline[0]=="@":
                input_format = "fastq"
            elif len(l)==1 and oneline[0]!=">" :
                input_format = "seq"
            elif len(l)==11:
                input_format = "qseq"
            elif oneline[0]==">" :
                input_format = "fasta"
            read_inf.close()


            #----------------------------------------------------------------
            read_id = ""
            seq = ""
            seq_ready = 0
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

                #---------------------------------------------------------------
                if seq_ready=="Y":
                    # Normalize the characters
                    seq=seq.upper().replace(".","N")

                    read_tag = [ m for m,n in [ (i, len(i)) for i in uniq(cut3_tag_lst)] if seq[0:n] == m ]
                    if len(read_tag) > 0 :
                        all_tagged += 1
                        for i in read_tag :
                            n_cut_tag_lst[i] += 1

                    seq = seq[(cut_s-1):cut_e] # cut_s start from 1 cycle by default

                    #-- Trimming adapter sequence ---

                    all_base_before_trim += len(seq)
                    if adapter != "" :
                        new_read = RemoveAdapter(seq, adapter, adapter_mismatch, Extra_base_cut_5end_adapter)
                        if len(new_read) < len(seq) :
                            all_tagged_trimmed += 1
                        seq = new_read
                    all_base_after_trim += len(seq)
                    if len(seq) <= 4 :
                        seq = "N" * (cut_e - cut_s)

                    # all reads will be considered, regardless of tags
                    #---------  trimmed_raw_BS_read and qscore ------------------
                    original_bs_reads[read_id] = seq
                    #---------  FW_C2T  ------------------
                    outf2.write('>%s\n%s\n'%(read_id, seq.replace('C', 'T')))
            fileinput.close()

            outf2.close()

            delete_files(read_file)
            logm("Processing input is done")
            #--------------------------------------------------------------------------------

            # mapping
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
            Multiple_hits=set()


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
                checking_genome_context = [output_genome[i:j] == k for i,j,k in zip(gx,gy,check_pattern) ]
                r_aln, g_aln = cigar_to_alignment(cigar, original_BS, origin_genome)

                if len(r_aln) == len(g_aln) :
                    my_region_serial, my_region_start, my_region_end = [-1, 0, 0]
                    if True in checking_genome_context :
                        try_pos = [mapped_location - len(i) for i,j in zip(cut5_tag_lst, checking_genome_context) if j][0]
                        my_region_serial, my_region_start, my_region_end = my_mappable_region(FW_regions[mapped_chr],
                                try_pos, FR)
                    if my_region_serial == 0 : # still be 0
                        # for some cases, read has no tags; searching the upstream sequence for tags
                        # print "[For debug]: FW read has no tags"
                        try_count = 0
                        try_pos = mapped_location - min_cut5_len + 1
                        while my_region_serial == 0 and try_count < MAX_TRY :
                            my_region_serial, my_region_start, my_region_end = my_mappable_region(FW_regions[mapped_chr],
                                                              try_pos, FR)
                            try_pos -= 1
                            try_count += 1

                        #if my_region_serial == 0 :
                        #    print "[For debug]: chr=", mapped_chr
                        #    print "[For debug]: +FW read still can not find fragment serial"
                        # Tip: sometimes "my_region_serial" is still 0 ...


                    N_mismatch = N_MIS(r_aln, g_aln)
#                    if N_mismatch <= int(max_mismatch_no) :
                    mm_no=float(max_mismatch_no)
                    if (mm_no>=1 and N_mismatch<=mm_no) or (mm_no<1 and N_mismatch<=(mm_no*len(r_aln)) ):
                        all_mapped_passed += 1
                        methy = methy_seq(r_aln, g_aln + next2bp)
                        mC_lst, uC_lst = mcounts(methy, mC_lst, uC_lst)
                        #---XS FILTER----------------
                        XS = 0
                        if XSteve :
                            if ('ZZZ' in methy.translate(None, "-XY")):
                                XS = 1
                            #
                        else :
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
                        all_base_mapped += len(original_BS)
                else :
                    print "[For debug]: reads not in same lengths"

            #print "start RC"
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
                #checking_genome_context = (output_genome[gx:gy] == check_pattern)
                checking_genome_context = [output_genome[i:j] == k for i,j,k in zip(gx,gy,check_pattern) ]
                r_aln, g_aln = cigar_to_alignment(cigar, original_BS, origin_genome)

                if len(r_aln) == len(g_aln) : # and checking_genome_context:
                    my_region_serial, my_region_start, my_region_end = [-1, 0, 0]
                    if True in checking_genome_context :
                        try_pos = [mapped_location + g_len - 1 + len(i) for i,j in zip(cut5_tag_lst, checking_genome_context) if j][0]
                        my_region_serial, my_region_start, my_region_end = my_mappable_region(RC_regions[mapped_chr],
                                                                           try_pos , FR)
                    if my_region_serial == 0 : # still be 0
                        # for some cases, read has no tags; searching the upstream sequence for tags
                        #print "[For debug]: RC Read has no tags"
                        try_count = 0
                        try_pos = mapped_location + g_len + min_cut5_len - 2
                        while my_region_serial == 0 and try_count < MAX_TRY :
                            my_region_serial, my_region_start, my_region_end = my_mappable_region(RC_regions[mapped_chr],
                                                                               try_pos, FR)
                            try_pos += 1
                            try_count += 1

                        #if my_region_serial == 0 :
                        #    print "[For debug]: chr=", mapped_chr
                        #    print "[For debug]: -FW read still cannot find fragment serial"


                    N_mismatch = N_MIS(r_aln, g_aln)
#                    if N_mismatch <= int(max_mismatch_no) :
                    mm_no=float(max_mismatch_no)
                    if (mm_no>=1 and N_mismatch<=mm_no) or (mm_no<1 and N_mismatch<=(mm_no*len(r_aln)) ):
                        all_mapped_passed += 1
                        methy = methy_seq(r_aln, g_aln + next2bp)
                        mC_lst, uC_lst = mcounts(methy, mC_lst, uC_lst)
                        #---XS FILTER----------------
                        XS = 0
                        if XSteve :
                            if ('ZZZ' in methy.translate(None, "-XY")):
                                XS = 1
                            #
                        else :
                            nCH = methy.count('y') + methy.count('z')
                            nmCH = methy.count('Y') + methy.count('Z')
                            if( (nmCH>XS_count) and nmCH/float(nCH+nmCH)>XS_pct ) :
                                XS = 1
                            #
                        #
                        num_mapped_RC_C2T += 1
                        outfile.store(header, N_mismatch, FR, mapped_chr, mapped_strand,
                                      mapped_location, cigar, original_BS, methy, XS,
                                      output_genome = output_genome,
                                      rrbs = True,
                                      my_region_serial = my_region_serial,
                                      my_region_start = my_region_start,
                                      my_region_end = my_region_end)
                        all_base_mapped += len(original_BS)
                else :
                    print "[For debug]: reads not in same lengths"


            # Finished both FW and RC
            logm("Done: %s (%d) \n" % (read_file, no_my_files))
            print "--> %s (%d) "%(read_file, no_my_files)
            del original_bs_reads
            delete_files(WC2T, CC2T)

    # End of directional library


    # ====================================================
    #  un-directional library
    # ====================================================

    elif asktag=="Y" :
        #----------------------------------------------------------------
        logm("Start reading and trimming the input sequences")

        input_fname = os.path.split(main_read_file)[1]
        for read_file in isplit_file(main_read_file, tmp_d(input_fname)+'-s-', no_small_lines):

            logm("Processing read file: %s" % read_file)
            original_bs_reads = {}
            no_my_files+=1
            random_id = ".tmp-"+str(random.randint(1000000,9999999))
            outfile2=tmp_d('Trimmed_C2T.fa'+random_id)
            outfile3=tmp_d('Trimmed_G2A.fa'+random_id)

            outf2=open(outfile2,'w')
            outf3=open(outfile3,'w')

            #--- Checking input format ------------------------------------------
            try :
                if read_file.endswith(".gz") : # support input file ending with ".gz"
                    read_inf = gzip.open(read_file, "rb")
                else :
                    read_inf=open(read_file,"r")
            except IOError:
                print "[Error] Cannot open input file : %s" % read_file
                exit(-1)

            oneline = read_inf.readline()
            if oneline == "" :
                oneline = "NNNN"
            l = oneline.split()
            input_format = ""
            if oneline[0]=="@":
                input_format = "fastq"
            elif len(l)==1 and oneline[0]!=">" :
                input_format = "seq"
            elif len(l)==11:
                input_format = "qseq"
            elif oneline[0]==">" :
                input_format = "fasta"
            read_inf.close()

            #----------------------------------------------------------------
            read_id = ""
            seq = ""
            seq_ready = 0
            line_no = 0
            for line in fileinput.input(read_file, openhook=fileinput.hook_compressed):
                if line == "" :
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

                #---------------------------------------------------------------
                if seq_ready=="Y":
                    # Normalize the characters
                    seq = seq.upper().replace(".","N")

                    read_tag = [ m for m,n in [ (i, len(i)) for i in uniq(cut3_tag_lst)] if seq[0:n] == m ]
                    if len(read_tag) > 0 :
                        all_tagged += 1
                        for i in read_tag :
                            n_cut_tag_lst[i] += 1

                    seq = seq[(cut_s-1):cut_e] # cut_s start from 1 cycle by default

                    #-- Trimming adapter sequence ---
                    all_base_before_trim += len(seq)
                    if adapter != "" :
                        new_read = RemoveAdapter(seq, adapter, adapter_mismatch, Extra_base_cut_5end_adapter)
                        if len(new_read) < len(seq) :
                            all_tagged_trimmed += 1
                        seq = new_read
                    all_base_after_trim += len(seq)
                    if len(seq) <= 4 :
                        seq = "N" * (cut_e - cut_s)

                    # all reads will be considered, regardless of tags
                    #---------  trimmed_raw_BS_read and qscore ------------------
                    original_bs_reads[read_id] = seq
                    #---------  FW_C2T  ------------------
                    outf2.write('>%s\n%s\n'%(read_id, seq.replace('C', 'T')))
                    #---------  RC_G2A  ------------------
                    outf3.write('>%s\n%s\n' % (read_id, seq.replace("G","A")))
            fileinput.close()

            outf2.close()

            delete_files(read_file)
            logm("Processing input is done")
            #--------------------------------------------------------------------------------

            # mapping
            #--------------------------------------------------------------------------------
            logm("Start mapping")
            WC2T = tmp_d("W_C2T_m"+str(max_mismatch_no)+".mapping"+random_id)
            CC2T = tmp_d("C_C2T_m"+str(max_mismatch_no)+".mapping"+random_id)
            WG2A = tmp_d("W_G2A_m"+str(max_mismatch_no)+".mapping"+random_id)
            CG2A = tmp_d("C_G2A_m"+str(max_mismatch_no)+".mapping"+random_id)

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

            FW_C2T_U,FW_C2T_R = extract_mapping(WC2T)
            RC_G2A_U,RC_G2A_R = extract_mapping(CG2A)

            FW_G2A_U,FW_G2A_R = extract_mapping(WG2A)
            RC_C2T_U,RC_C2T_R = extract_mapping(CC2T)

            logm("Extracting alignments is done")

            #----------------------------------------------------------------
            # get unique-hit reads
            #----------------------------------------------------------------
            Union_set = set(FW_C2T_U.iterkeys()) | set(RC_G2A_U.iterkeys()) | set(FW_G2A_U.iterkeys()) | set(RC_C2T_U.iterkeys())

            Unique_FW_C2T = set() # +
            Unique_RC_G2A = set() # +
            Unique_FW_G2A = set() # -
            Unique_RC_C2T = set() # -
            Multiple_hits = set()

            for x in Union_set :
                _list = []
                for dx in [FW_C2T_U, RC_G2A_U, FW_G2A_U, RC_C2T_U]:
                    mis_lst = dx.get(x,[99])
                    mis = int(mis_lst[0])
                    _list.append(mis)
                for dx in [FW_C2T_R, RC_G2A_R, FW_G2A_R, RC_C2T_R]:
                    mis = dx.get(x,99)
                    _list.append(mis)
                mini = min(_list)
                if _list.count(mini) == 1:
                    mini_index = _list.index(mini)
                    if mini_index == 0:
                        Unique_FW_C2T.add(x)
                    elif mini_index == 1:
                        Unique_RC_G2A.add(x)
                    elif mini_index == 2:
                        Unique_FW_G2A.add(x)
                    elif mini_index == 3:
                        Unique_RC_C2T.add(x)
                    else :
                        Multiple_hits.add(x)
                else :
                    Multiple_hits.add(x)
            # write reads rejected by Multiple Hits to file
            if show_multiple_hit is not None :
                #outf_MH = open(show_multiple_hit,'w')
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

            FW_C2T_uniq_lst = [[FW_C2T_U[u][1],u] for u in Unique_FW_C2T]
            FW_G2A_uniq_lst = [[FW_G2A_U[u][1],u] for u in Unique_FW_G2A]
            RC_C2T_uniq_lst = [[RC_C2T_U[u][1],u] for u in Unique_RC_C2T]
            RC_G2A_uniq_lst = [[RC_G2A_U[u][1],u] for u in Unique_RC_G2A]
            FW_C2T_uniq_lst.sort()
            RC_C2T_uniq_lst.sort()
            FW_G2A_uniq_lst.sort()
            RC_G2A_uniq_lst.sort()
            FW_C2T_uniq_lst = [x[1] for x in FW_C2T_uniq_lst]
            RC_C2T_uniq_lst = [x[1] for x in RC_C2T_uniq_lst]
            FW_G2A_uniq_lst = [x[1] for x in FW_G2A_uniq_lst]
            RC_G2A_uniq_lst = [x[1] for x in RC_G2A_uniq_lst]

            del Unique_FW_C2T
            del Unique_FW_G2A
            del Unique_RC_C2T
            del Unique_RC_G2A


            #----------------------------------------------------------------
            # Post-filtering reads
            # ---- FW_C2T  ---- un-directional
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
                checking_genome_context = [output_genome[i:j] == k for i,j,k in zip(gx,gy,check_pattern) ]
                r_aln, g_aln = cigar_to_alignment(cigar, original_BS, origin_genome)

                if len(r_aln) == len(g_aln) :
                    my_region_serial, my_region_start, my_region_end = [-1, 0, 0]
                    if True in checking_genome_context :
                        try_pos = [mapped_location - len(i) for i,j in zip(cut5_tag_lst, checking_genome_context) if j][0]
                        my_region_serial, my_region_start, my_region_end = my_mappable_region(FW_regions[mapped_chr],
                                try_pos, FR)
                    if my_region_serial == 0 : # still be 0
                        # for some cases, read has no tags; searching the upstream sequence for tags
                        # print "[For debug]: FW read has no tags"
                        try_count = 0
                        try_pos = mapped_location - min_cut5_len + 1
                        while my_region_serial == 0 and try_count < MAX_TRY :
                            my_region_serial, my_region_start, my_region_end = my_mappable_region(FW_regions[mapped_chr],
                                                              try_pos, FR)
                            try_pos -= 1
                            try_count += 1

                    N_mismatch = N_MIS(r_aln, g_aln)
#                    if N_mismatch <= int(max_mismatch_no) :
                    mm_no=float(max_mismatch_no)
                    if (mm_no>=1 and N_mismatch<=mm_no) or (mm_no<1 and N_mismatch<=(mm_no*len(r_aln)) ):
                        all_mapped_passed += 1
                        methy = methy_seq(r_aln, g_aln + next2bp)
                        mC_lst, uC_lst = mcounts(methy, mC_lst, uC_lst)
                        # ---XS FILTER----------------
                        XS = 0
                        if XSteve:
                            if ('ZZZ' in methy.translate(None, "-XY")):
                                XS = 1
                                #
                        else:
                            nCH = methy.count('y') + methy.count('z')
                            nmCH = methy.count('Y') + methy.count('Z')
                            if ((nmCH > XS_count) and nmCH / float(nCH + nmCH) > XS_pct):
                                XS = 1
                                #
                        #
                        num_mapped_FW_C2T += 1
                        outfile.store(header, N_mismatch, FR, mapped_chr, mapped_strand,
                                      mapped_location, cigar, original_BS, methy, XS,
                                      output_genome = output_genome,
                                      rrbs = True,
                                      my_region_serial = my_region_serial,
                                      my_region_start = my_region_start,
                                      my_region_end = my_region_end)
                        all_base_mapped += len(original_BS)
                else :
                    print "[For debug]: reads not in same lengths"


            # ---- RC_C2T ---- un-directional
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
                checking_genome_context = [output_genome[i:j] == k for i,j,k in zip(gx,gy,check_pattern) ]
                r_aln, g_aln = cigar_to_alignment(cigar, original_BS, origin_genome)

                if len(r_aln) == len(g_aln) : # and checking_genome_context:
                    my_region_serial, my_region_start, my_region_end = [-1, 0, 0]
                    if True in checking_genome_context :
                        try_pos = [mapped_location + g_len - 1 + len(i) for i,j in zip(cut5_tag_lst, checking_genome_context) if j][0]
                        my_region_serial, my_region_start, my_region_end = my_mappable_region(RC_regions[mapped_chr],
                                                                           try_pos , FR)
                    if my_region_serial == 0 : # still be 0
                        # for some cases, read has no tags; searching the upstream sequence for tags
                        #print "[For debug]: RC Read has no tags"
                        try_count = 0
                        try_pos = mapped_location + g_len + min_cut5_len - 2
                        while my_region_serial == 0 and try_count < MAX_TRY :
                            my_region_serial, my_region_start, my_region_end = my_mappable_region(RC_regions[mapped_chr],
                                                                               try_pos, FR)
                            try_pos += 1
                            try_count += 1

                    N_mismatch = N_MIS(r_aln, g_aln)
#                    if N_mismatch <= int(max_mismatch_no) :
                    mm_no=float(max_mismatch_no)
                    if (mm_no>=1 and N_mismatch<=mm_no) or (mm_no<1 and N_mismatch<=(mm_no*len(r_aln)) ):
                        all_mapped_passed += 1
                        methy = methy_seq(r_aln, g_aln + next2bp)
                        mC_lst, uC_lst = mcounts(methy, mC_lst, uC_lst)
                        #---XS FILTER----------------
                        XS = 0
                        if XSteve:
                            if ('ZZZ' in methy.translate(None, "-XY")):
                                XS = 1
                                #
                        else:
                            nCH = methy.count('y') + methy.count('z')
                            nmCH = methy.count('Y') + methy.count('Z')
                            if ((nmCH > XS_count) and nmCH / float(nCH + nmCH) > XS_pct):
                                XS = 1
                                #
                        #
                        num_mapped_RC_C2T += 1
                        outfile.store(header, N_mismatch, FR, mapped_chr, mapped_strand,
                                      mapped_location, cigar, original_BS, methy, XS,
                                      output_genome = output_genome,
                                      rrbs = True,
                                      my_region_serial = my_region_serial,
                                      my_region_start = my_region_start,
                                      my_region_end = my_region_end)
                        all_base_mapped += len(original_BS)

                else :
                    print "[For debug]: reads not in same lengths"


            # ---- FW_G2A  ---- un-directional
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
                checking_genome_context = [output_genome[i:j] == k for i,j,k in zip(gx,gy,check_pattern) ]
                r_aln, g_aln = cigar_to_alignment(cigar, original_BS, origin_genome)

                if len(r_aln) == len(g_aln) :
                    my_region_serial, my_region_start, my_region_end = [-1, 0, 0]
                    if True in checking_genome_context :
                        try_pos = [mapped_location - len(i) for i,j in zip(cut5_tag_lst, checking_genome_context) if j][0]
                        my_region_serial, my_region_start, my_region_end = my_mappable_region(FW_regions[mapped_chr],
                                                            try_pos, FR)
                    if my_region_serial == 0 : # still be 0
                        # for some cases, read has no tags; searching the upstream sequence for tags
                        #print "[For debug]: FW read has no tags"
                        try_count = 0
                        try_pos = mapped_location - min_cut5_len + 1
                        while my_region_serial == 0 and try_count < MAX_TRY :
                            my_region_serial, my_region_start, my_region_end = my_mappable_region(FW_regions[mapped_chr],
                                            try_pos, FR)
                            try_pos += 1
                            try_count += 1
                        #if my_region_serial == 0 :
                        #    print "[For debug]: chr=", mapped_chr
                        #    print "[For debug]: FW_G2A read still can not find fragment serial"
                        # Tip: sometimes "my_region_serial" is still 0 ...

                    N_mismatch = N_MIS(r_aln, g_aln)
#                    if N_mismatch <= int(max_mismatch_no) :
                    max_mismatch_no=float(max_mismatch_no)
                    if (mm_no>=1 and N_mismatch<=mm_no) or (mm_no<1 and N_mismatch<=(mm_no*len(r_aln)) ):
                        all_mapped_passed += 1
                        methy = methy_seq(r_aln, g_aln + next2bp)
                        mC_lst, uC_lst = mcounts(methy, mC_lst, uC_lst)
                        #---XS FILTER----------------
                        XS = 0
                        # ---XS FILTER----------------
                        XS = 0
                        if XSteve:
                            if ('ZZZ' in methy.translate(None, "-XY")):
                                XS = 1
                                #
                        else:
                            nCH = methy.count('y') + methy.count('z')
                            nmCH = methy.count('Y') + methy.count('Z')
                            if ((nmCH > XS_count) and nmCH / float(nCH + nmCH) > XS_pct):
                                XS = 1
                                #
                        #
                        num_mapped_FW_G2A += 1
                        outfile.store(header, N_mismatch, FR, mapped_chr, mapped_strand,
                                      mapped_location, cigar, original_BS, methy, XS,
                                      output_genome = output_genome,
                                      rrbs = True,
                                      my_region_serial = my_region_serial,
                                      my_region_start = my_region_start,
                                      my_region_end = my_region_end)
                        all_base_mapped += len(original_BS)
                else :
                    print "[For debug]: reads not in same lengths"


            # ---- RC_G2A ---- un-directional
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
                checking_genome_context = [output_genome[i:j] == k for i,j,k in zip(gx,gy,check_pattern) ]
                r_aln, g_aln = cigar_to_alignment(cigar, original_BS, origin_genome)

                if len(r_aln) == len(g_aln) : # and checking_genome_context:
                    my_region_serial, my_region_start, my_region_end = [-1, 0, 0]
                    if True in checking_genome_context :
                        try_pos = [mapped_location + g_len - 1 + len(i) for i,j in zip(cut5_tag_lst, checking_genome_context) if j][0]
                        my_region_serial, my_region_start, my_region_end = my_mappable_region(RC_regions[mapped_chr],
                                                                    mapped_location + g_len + min_cut5_len -1, FR)
                    if my_region_serial == 0 : # still be 0
                        # for some cases, read has no tags; searching the upstream sequence for tags
                        #print "[For debug]: RC Read has no tags"
                        try_count = 0
                        try_pos = mapped_location + g_len + min_cut5_len - 2
                        while try_count < MAX_TRY :
                            my_region_serial, my_region_start, my_region_end = my_mappable_region(RC_regions[mapped_chr],
                                                    try_pos, FR)
                            try_pos += 1
                            try_count += 1
                        #
                        #if my_region_serial == 0 :
                        #    print "[For debug]: chr=", mapped_chr
                        #    print "[For debug]: RC_C2A read still cannot find fragment serial"
                    #
                    N_mismatch = N_MIS(r_aln, g_aln)
                    #if N_mismatch <= int(max_mismatch_no) :
                    mm_no=float(max_mismatch_no)
                    if (mm_no>=1 and N_mismatch<=mm_no) or (mm_no<1 and N_mismatch<=(mm_no*len(r_aln)) ):
                        all_mapped_passed += 1
                        methy = methy_seq(r_aln, g_aln + next2bp)
                        mC_lst, uC_lst = mcounts(methy, mC_lst, uC_lst)
                        # ---XS FILTER----------------
                        XS = 0
                        if XSteve:
                            if ('ZZZ' in methy.translate(None, "-XY")):
                                XS = 1
                                #
                        else:
                            nCH = methy.count('y') + methy.count('z')
                            nmCH = methy.count('Y') + methy.count('Z')
                            if ((nmCH > XS_count) and nmCH / float(nCH + nmCH) > XS_pct):
                                XS = 1
                                #
                        #
                        num_mapped_RC_G2A += 1
                        outfile.store(header, N_mismatch, FR, mapped_chr, mapped_strand,
                                      mapped_location, cigar, original_BS, methy, XS,
                                      output_genome = output_genome,
                                      rrbs = True,
                                      my_region_serial = my_region_serial,
                                      my_region_start = my_region_start,
                                      my_region_end = my_region_end)
                        all_base_mapped += len(original_BS)
                else :
                    print "[For debug]: reads not in same lengths"

            # Finished both FW and RC
            logm("Done: %s (%d) \n" % (read_file, no_my_files))
            print "--> %s (%d) "%(read_file, no_my_files)
            del original_bs_reads
            delete_files(WC2T, CC2T, WG2A, CG2A)
    # End of un-directional library

    delete_files(tmp_path)

    if show_multiple_hit is not None :
        outf_MH.close()
    if show_unmapped_hit is not None :
        outf_UH.close()

    logm("----------------------------------------------")
    logm("Number of raw reads: %d " % all_raw_reads)
    if all_raw_reads>0:
        tag_percent = (float(all_tagged)/all_raw_reads) if all_raw_reads>0 else 0
        logm("Number of raw reads with CGG/TGG at 5' end: %d (%1.3f)" % (all_tagged, tag_percent) )
        for kk in range(len(n_cut_tag_lst)):
            tag_percent = (float(n_cut_tag_lst[cut3_tag_lst[kk]])/all_raw_reads) if all_raw_reads>0 else 0
            logm("Number of raw reads with tag %s: %d (%1.3f)" % (cut3_tag_lst[kk],n_cut_tag_lst[cut3_tag_lst[kk]], tag_percent) )
        logm("Number of bases in total: %d " % all_base_before_trim)
        if adapter!="" :
            logm("Number of reads having adapter removed: %d " % all_tagged_trimmed)
            trim_percent = (float(all_base_after_trim)/all_base_before_trim ) if all_base_before_trim>0 else 0
            logm("Number of bases after trimming the adapters: %d (%1.3f)" % (all_base_after_trim, trim_percent ) )
        logm("Number of reads are rejected because of multiple hits: %d\n" % len(Multiple_hits) )
        logm("Number of unique-hits reads (before post-filtering): %d" % all_mapped)
        logm("  %d uniquely aligned reads, passed fragment check, with mismatches <= %s"%(all_mapped_passed, max_mismatch_no))
        Mappability = (100*float(all_mapped_passed)/all_raw_reads) if all_raw_reads>0 else 0
        logm("Mappability = %1.4f%%" % Mappability )
        #logm("Total bases of uniquely mapped reads %7d"% all_base_mapped )
        if asktag == "Y": # un-directional
            logm("  %7d FW reads mapped to Watson strand" % (num_mapped_FW_C2T) )
            logm("  %7d RC reads mapped to Watson strand" % (num_mapped_FW_G2A) )
            logm("  %7d FW reads mapped to Crick strand" % (num_mapped_RC_C2T) )
            logm("  %7d RC reads mapped to Crick strand" % (num_mapped_RC_G2A) )
            # the variable name 'num_mapped_RC_G2A' seems not consistent with illustration
            # according to literal meaning
        elif asktag == "N": # directional
            logm("  %7d FW reads mapped to Watson strand" % (num_mapped_FW_C2T) )
            logm("  %7d FW reads mapped to Crick strand" % (num_mapped_RC_C2T) )
        logm("Total bases of uniquely mapped reads : %7d" % all_base_mapped )
        #
        n_CG  = mC_lst[0] + uC_lst[0]
        n_CHG = mC_lst[1] + uC_lst[1]
        n_CHH = mC_lst[2] + uC_lst[2]
        #
        logm("----------------------------------------------")
        logm("Methylated C in mapped reads ")
        logm(" mCG  %1.3f%%" % ( (100 * float(mC_lst[0]) / n_CG) if n_CG!=0 else 0))
        logm(" mCHG %1.3f%%" % ( (100 * float(mC_lst[1]) / n_CHG) if n_CHG!=0 else 0))
        logm(" mCHH %1.3f%%" % ( (100 * float(mC_lst[2]) / n_CHH) if n_CHH!=0 else 0))
        #
    logm("-------------------------------- " )
    logm("File : %s" % main_read_file )
    elapsed("Resource / CPU time")
    logm("------------------- END ----------------------")
    close_log()


