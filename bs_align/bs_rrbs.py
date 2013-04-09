import fileinput, random, math, os.path
from bs_index.rrbs_build import FWD_MAPPABLE_REGIONS, REV_MAPPABLE_REGIONS
from bs_utils.utils import *

from bs_align.bs_single_end import extract_mapping
from bs_align_utils import *

def my_mapable_region(chr_regions, mapped_location, FR): # start_position (first C), end_position (last G), serial, sequence
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
    elif FR == "-FW":
        my_location = str(mapped_location-1)
        if my_location in chr_regions:
            my_lst = chr_regions[my_location]
            out_end = int(my_location)
            out_start = my_lst[0]
            out_serial = my_lst[1]
    return out_serial, out_start, out_end



#----------------------------------------------------------------

def bs_rrbs(main_read_file, mytag, adapter_file, cut1, cut2, no_small_lines, indexname, aligner_command, db_path, tmp_path, outfile, XS_pct, XS_count, adapter_mismatch):
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
        except IOError:
            print "[Error] Cannot find adapter file : %s !" % adapter_file
            exit(-1)
        adapter_seq = adapter_inf.readline().strip()
        adapter_inf.close()

    #----------------------------------------------------------------
    print "Read filename: %s" % main_read_file
    print "Starting Msp-1 tag: %s"% mytag
    print "The last cycle (for mapping): %d"% cut2
    print "Bowtie path: %s" % aligner_command
    print "Reference genome library path: %s" % db_path
    print "Number of mismatches allowed: %s" % indexname

    logm("I Read filename: %s" % main_read_file)
    logm("I Starting Msp-1 tag: %s" % mytag )
    logm("I The last cycle (for mapping): %d" % cut2 )
    logm("I Bowtie path: %s" % aligner_command )
    logm("I Reference genome library path: %s" % db_path )
    logm("I Number of mismatches allowed: %s" % indexname)
    logm("I adapter seq: %s" % adapter_seq)
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
        input_format=""
        if oneline[0]=="@":	# Illumina GAII FastQ (Lister et al Nature 2009)
            input_format="FastQ"
            n_fastq=0
        elif len(l)==1 and oneline[0]!=">": 	# pure sequences
            input_format="list of sequences"
        elif len(l)==11:	# Illumina GAII qseq file
            input_format="Illumina GAII qseq file"
        elif oneline[0]==">":	# fasta
            input_format="fasta"
            n_fasta=0
        read_inf.close()

        #----------------------------------------------------------------
        id=""
        seq=""
        seq_ready=0
        for line in fileinput.input(read_file):
            l=line.split()

            if input_format=="Old Solexa Seq file":
                all_raw_reads+=1
                id=str(all_raw_reads)
                id=id.zfill(12)
                seq=l[4]
                seq_ready="Y"
            elif input_format=="list of sequences":
                all_raw_reads+=1
                id=str(all_raw_reads)
                id=id.zfill(12)
                seq=l[0]
                seq_ready="Y"
            elif input_format=="FastQ":
                m_fastq=math.fmod(n_fastq,4)
                n_fastq+=1
                seq_ready="N"
                if m_fastq==0:
                    all_raw_reads+=1
                    id=str(all_raw_reads)
                    id=id.zfill(12)
                    seq=""
                elif m_fastq==1:
                    seq=l[0]
                    seq_ready="Y"
                else:
                    seq=""
            elif input_format=="Illumina GAII qseq file":
                all_raw_reads+=1
                id=str(all_raw_reads)
                id=id.zfill(12)
                seq=l[8]
                seq_ready="Y"
            elif input_format=="fasta":
                m_fasta=math.fmod(n_fasta,2)
                n_fasta+=1
                seq_ready="N"
                if m_fasta==0:
                    all_raw_reads+=1
                    id=l[0][1:]
                    seq=""
                elif m_fasta==1:
                    seq=l[0]
                    seq_ready="Y"
                else:
                    seq=""
            #---------------------------------------------------------------
            if seq_ready=="Y":
                seq=seq[0:cut2].upper().replace(".","N") #<----------------------selecting 0..52 from 1..72  -e 52

                #-- Selecting Reads with mytag (i.e., CGG or TGG or CGA) -----------------------
                has_tag="N"
                for i in range(cut2):
                    if seq[i:i+3] in mytag_lst and has_tag=="N":
                        all_tagged+=1
                        n_mytag_lst[seq[i:i+3]]+=1
                        has_tag="Y"
                        seq=seq[i:]

                        #-- Trimming adapter sequence ---
                        # New way to remove the adaptor
                        if adapter_seq != "":
                            new_read = RemoveAdapter(seq, adapter_seq, adapter_mismatch)
                           # new_read = RemoveAdapter(seq, adapter_seq, 1)
                            if len(new_read) < len(seq) :
                                all_tagged_trimed += 1
                            seq = new_read
                        if len(seq)<=4:
                            seq = "N" * cut2

                        break
                    #print "seq=", seq
                if has_tag=="Y":
                    #---------  trimmed_raw_BS_read and qscore ------------------
                    original_bs_reads[id] = seq

                    #---------  FW_C2T  ------------------
                    outf2.write('>%s\n%s\n'%(id, seq.replace('C', 'T')))
                #print seq
        # for_end
        fileinput.close()

        outf2.close()

        delete_files(read_file)
        logm("Processing input is done")
        #--------------------------------------------------------------------------------
        # Bowtie mapping
        #--------------------------------------------------------------------------------
        WC2T=tmp_d("W_C2T_m"+indexname+".mapping"+random_id)
        CC2T=tmp_d("C_C2T_m"+indexname+".mapping"+random_id)

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
            list=[]
            for dx in [FW_C2T_U,RC_C2T_U]:
                mis_lst=dx.get(x,[99])
                mis=int(mis_lst[0])
                list.append(mis)
            for dx in [FW_C2T_R,RC_C2T_R]:
                mis=dx.get(x,99)
                list.append(mis)
            mini=min(list)
            if list.count(mini)==1:
                mini_index=list.index(mini)
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

        nn=0
        for ali_unique_lst, ali_dic in [(FW_uniq_lst, FW_C2T_U), (RC_uniq_lst, RC_C2T_U)]:
            nn += 1
            mapped_chr0 = ""
            for header in ali_unique_lst:
                _, mapped_chr, mapped_location, cigar = ali_dic[header]
                original_BS = original_bs_reads[header]
                #-------------------------------------
                if mapped_chr != mapped_chr0:
                    FW_chr_regions = deserialize(db_d(FWD_MAPPABLE_REGIONS(mapped_chr)))
                    RC_chr_regions = deserialize(db_d(REV_MAPPABLE_REGIONS(mapped_chr)))
                    my_gseq = deserialize(db_d(mapped_chr))
                    chr_length = len(my_gseq)
                    mapped_chr0=mapped_chr
                r_start, r_end, g_len = get_read_start_end_and_genome_length(cigar)
                all_mapped+=1

#                checking_first_C = False
                if nn == 1: 							# +FW mapped to + strand:
                    FR = "+FW"
#                    mapped_location += 1 # 1 based (after plus 1)
#                    origin_genome_long = my_gseq[mapped_location - 2 - 1 : mapped_location + g_len + 2 - 1]
#                    checking_first_C = (origin_genome_long[1:5] == "CCGG")
                    mapped_strand = "+"
#                    origin_genome = origin_genome_long[2:-2]

                elif nn==2: 						# -FW mapped to - strand:
                    mapped_strand = "-"
                    FR = "-FW"
                    mapped_location = chr_length - mapped_location - g_len
#                    origin_genome_long = my_gseq[mapped_location - 2 - 1 : mapped_location + g_len + 2 - 1]
#                    origin_genome_long = reverse_compl_seq(origin_genome_long)
#                    checking_first_C = (origin_genome_long[1:5] == "CCGG")
#                    origin_genome = origin_genome_long[2:-2]


                origin_genome, next, output_genome = get_genomic_sequence(my_gseq, mapped_location, mapped_location + g_len, mapped_strand)
                checking_first_C = (output_genome[1:6] == "C_CGG")

                r_aln, g_aln = cigar_to_alignment(cigar, original_BS, origin_genome)

                if len(r_aln) == len(g_aln) and checking_first_C:
                    #---------------------------------------------
                    if FR=="+FW":
                        my_region_serial, my_region_start, my_region_end = my_mapable_region(FW_chr_regions, mapped_location + 1, "+FW")
                    elif FR=="-FW":
                        my_region_serial, my_region_start, my_region_end = my_mapable_region(RC_chr_regions, mapped_location + g_len + 1, "-FW")
                    #---------------------------------------------

                    N_mismatch = N_MIS(r_aln, g_aln) #+ original_BS_length - (r_end - r_start) # mismatches in the alignment + soft clipped nucleotides

                    if N_mismatch <= int(indexname) and my_region_serial != 0:
                        all_mapped_passed += 1
                        #---------------------------------------------

#                        output_genome = "%s_%s_%s" % (origin_genome_long[0:2], origin_genome, origin_genome_long[-2:])

                        methy = methy_seq(r_aln, g_aln + next)

                        mC_lst, uC_lst = mcounts(methy, mC_lst, uC_lst)

                        #---XS FILTER----------------
                        #XS = 1 if "ZZZ" in methy.replace('-', '') else 0
                        XS = 0
                        nCH = methy.count('y') + methy.count('z')
                        nmCH = methy.count('Y') + methy.count('Z')
                        # print nCH, nmCH
                        if( (nmCH>XS_count) and nmCH/float(nCH+nmCH)>XS_pct ) :
                            # print "find one", XS
                            XS = 1

                        outfile.store(header, N_mismatch, FR, mapped_chr, mapped_strand, mapped_location, cigar, original_BS, methy, XS, output_genome = output_genome, rrbs = True, my_region_serial = my_region_serial, my_region_start = my_region_start, my_region_end = my_region_end)


#        logm("Done: %s (%d/%d) \n" % (read_file, no_my_files, len(my_files)))
#        print "--> %s (%d/%d) "%(read_file, no_my_files, len(my_files))
        logm("Done: %s (%d) \n" % (read_file, no_my_files))
        print "--> %s (%d) "%(read_file, no_my_files)
        del original_bs_reads
        delete_files(WC2T, CC2T)

    delete_files(tmp_path)

    #----------------------------------------------------------------


    logm("O Number of raw reads: %d "% all_raw_reads)
    if all_raw_reads >0:
        logm("O Number of CGG/TGG tagged reads: %d (%1.3f)"%(all_tagged,float(all_tagged)/all_raw_reads))
        for kk in range(len(n_mytag_lst)):
            logm("O Number of raw reads with %s tag: %d (%1.3f)"%(mytag_lst[kk],n_mytag_lst[mytag_lst[kk]],float(n_mytag_lst[mytag_lst[kk]])/all_raw_reads))
        logm("O Number of CGG/TGG reads having adapter removed: %d "%all_tagged_trimed)
        logm("O Number of unique-hits reads for post-filtering: %d"%all_mapped)

        logm("O ------ %d uniquely aligned reads, passed fragment check, with mismatches <= %s"%(all_mapped_passed, indexname))
        logm("O Mapability= %1.4f%%"%(100*float(all_mapped_passed)/all_raw_reads))

        n_CG=mC_lst[0]+uC_lst[0]
        n_CHG=mC_lst[1]+uC_lst[1]
        n_CHH=mC_lst[2]+uC_lst[2]

        logm("----------------------------------------------")
        logm("M Methylated C in mapped reads ")
        logm("M mCG %1.3f%%"%(100*float(mC_lst[0])/n_CG))
        logm("M mCHG %1.3f%%"%(100*float(mC_lst[1])/n_CHG))
        logm("M mCHH %1.3f%%"%(100*float(mC_lst[2])/n_CHH))

    logm("----------------------------------------------")
    logm("------------------- END ----------------------")
    elapsed(main_read_file)

    close_log()
