BS-Seeker2
=========

BS-Seeker2 (BS Seeker 2) performs accurate and fast mapping of bisulfite-treated short reads. BS-Seeker2 is an updated version on BS-Seeker.

0. Availability
============

Homepage of [BS-Seeker2](http://pellegrini.mcdb.ucla.edu/BS_Seeker2/).

The source code for this package is available from
[https://github.com/BSSeeker/BSseeker2](https://github.com/BSSeeker/BSseeker2).
Also, you can use an instance of BS-Seeker 2 in Galaxy from [http://galaxy.hoffman2.idre.ucla.edu](http://galaxy.hoffman2.idre.ucla.edu).
(Label: "NGS: Methylation Mapping"/"Methylation Map with BS Seeker2")


1. Remarkable new features
============
* Reduced index for RRBS, accelerating the mapping speed and increasing mappability
* Allowing local alignment with Bowtie 2, increased the mappability

2. Other features
============
* Supported library types
	- whole genome-wide bisulfite sequencing (WGBS)
	- reduced representative bisulfite sequencing (RRBS)

* Supported formats for input file
	- [fasta](http://en.wikipedia.org/wiki/FASTA_format)
	- [fastq](http://en.wikipedia.org/wiki/FASTQ_format)
	- [qseq](http://jumpgate.caltech.edu/wiki/QSeq)
	- pure sequence (one-line one-sequence)

* Supported alignment tools
	- [bowtie](http://bowtie-bio.sourceforge.net/index.shtml) : Single-seed
	- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) : Multiple-seed, gapped-alignment
		- [local alignment](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#local-alignment-example)
		- [end-to-end alignment](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#end-to-end-alignment-example)
	- [SOAP](http://soap.genomics.org.cn/)

* Supported formats for mapping results
	- [BAM](http://genome.ucsc.edu/FAQ/FAQformat.html#format5.1)
	- [SAM](http://samtools.sourceforge.net/)
	- [BS-seeker 1](http://pellegrini.mcdb.ucla.edu/BS_Seeker/USAGE.html)

3. System requirements
============

* Linux or Mac OS platform
* One of the following Aligner
  - [bowtie](http://bowtie-bio.sourceforge.net/) (fast)
  - [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/) (Default)
  - [soap](http://soap.genomics.org.cn/)
  - [rmap](http://www.cmb.usc.edu/people/andrewds/rmap/)
* [Python](http://www.python.org/download/) (Version 2.6 +)

  (It is normally pre-installed in Linux. Type " python -V" to see the installed version.)

* [pysam](http://code.google.com/p/pysam/) package is needed.

  (Read "Questions & Answers" if you have problem when installing this package.)



4. Modules' descriptions
============

(0) FilterReads.py
------------

Optional and independent module.
Some reads would be extremely amplified during the PCR. This script helps you get unique reads before doing the mapping. You can decide whether or not to filter reads before doing the mapping.

##Usage :

	$ python FilterReads.py
	Usage: FilterReads.py -i <input> -o <output> [-k]
	Author : Guo, Weilong; 2012-11-10
	Unique reads for qseq/fastq/fasta/sequencce, and filter
	low quality file in qseq file.

	Options:
	  -h, --help  show this help message and exit
	  -i FILE     Name of the input qseq/fastq/fasta/sequence file
	  -o FILE     Name of the output file
	  -k          Would not filter low quality reads if specified


##Tip :

- This step is not suggested for RRBS library, as reads from RRBS library would more likely from the same location.


(1) bs_seeker2-build.py
------------

Module to build the index for BS-Seeker2.


##Usage :

    $ python bs_seeker2-build.py -h
    Usage: bs_seeker2-build.py [options]

    Options:
      -h, --help            show this help message and exit
      -f FILE, --file=FILE  Input your reference genome file (fasta)
      --aligner=ALIGNER     Aligner program to perform the analysis: bowtie,
                            bowtie2, soap, rmap [Default: bowtie]
      -p PATH, --path=PATH  Path to the aligner program. Detected:
                            bowtie: ~/install/bowtie
                            bowtie2: ~/install/bowtie2
                            rmap: ~/install/rmap_/bin
                            soap: ~/install/soap/
      -d DBPATH, --db=DBPATH
                            Path to the reference genome library (generated in
                            preprocessing genome) [Default: ~/install
                            /BSseeker2/bs_utils/reference_genomes]
      -v, --version         show version of BS-Seeker2

      Reduced Representation Bisulfite Sequencing Options:
        Use this options with conjuction of -r [--rrbs]

        -r, --rrbs          Build index specially for Reduced Representation
                            Bisulfite Sequencing experiments. Genome other than
                            certain fragments will be masked. [Default: False]
        -l LOW_BOUND, --low=LOW_BOUND
                            lower bound of fragment length (excluding recognition
                            sequence such as C-CGG) [Default: 20]
        -u UP_BOUND, --up=UP_BOUND
                            upper bound of fragment length (excluding recognition
                            sequence such as C-CGG ends) [Default: 500]
        -c CUT_FORMAT, --cut-site=CUT_FORMAT
                            Cut sites of restriction enzyme. Ex: MspI(C-CGG),
                            Mael:(C-TAG), double-enzyme MspI&Mael:(C-CGG,C-TAG).
                            [Default: C-CGG]


##Example

* Build genome index for WGBS using bowtie, path of bowtie should be included in $PATH

        python bs_seeker2-build.py -f genome.fa --aligner=bowtie

* Build genome index for RRBS with default parameters specifying the path for bowtie2

        python bs_seeker2-build.py -f genome.fa --aligner=bowtie2 -p ~/install/bowtie2-2.0.0-beta7/ -r

* Build genome index for RRBS library using bowite2, with fragment lengths ranging [40bp, 400bp]

        python bs_seeker2-build.py -f genome.fa -r -l 40 -u 400 --aligner=bowtie2

* Build genome index for RRBS library for double-enzyme : MspI (C-CGG) & ApeKI (G-CWGC, where W=A|T, see [IUPAC code](http://www.bioinformatics.org/sms/iupac.html))

        python bs_seeker2-build.py -f genome.fa -r -c C-CGG,G-CWGC --aligner=bowtie

##Tips:

- Index built for BS-Seeker2 is different from the index for BS-Seeker 1.
For RRBS, you need to specify "-r" in the parameters. Also, you need to specify LOW_BOUND and UP_BOUND for the range of fragment lengths according your protocol.

- The fragment length is different from read length. Fragments refers to the DNA fragments which you get by size-selection step (i.e. gel-cut oor AMPure beads). Lengths of fragments are supposed to be in a range, such as [50bp,250bp].

- The indexes for RRBS and WGBS are different. Also, indexes for RRBS are specific for fragment length parameters (LOW_BOUND and UP_BOUND).




(2) bs_seeker2-align.py
------------

Module to map reads on 3-letter converted genome.

##Usage :

    $ python ~/install/BSseeker2/bs_seeker2-align.py -h
    Usage: bs_seeker2-align.py {-i <single> | -1 <mate1> -2 <mate2>} -g <genome.fa> [options]

    Options:
      -h, --help            show this help message and exit

      For single end reads:
        -i INFILE, --input=INFILE
                            Input read file (FORMAT: sequences, qseq, fasta,
                            fastq). Ex: read.fa or read.fa.gz

      For pair end reads:
        -1 FILE, --input_1=FILE
                            Input read file, mate 1 (FORMAT: sequences, qseq,
                            fasta, fastq)
        -2 FILE, --input_2=FILE
                            Input read file, mate 2 (FORMAT: sequences, qseq,
                            fasta, fastq)
        -I MIN_INSERT_SIZE, --minins=MIN_INSERT_SIZE
                            The minimum insert size for valid paired-end
                            alignments [Default: 0]
        -X MAX_INSERT_SIZE, --maxins=MAX_INSERT_SIZE
                            The maximum insert size for valid paired-end
                            alignments [Default: 500]

      Reduced Representation Bisulfite Sequencing Options:
        -r, --rrbs          Map reads to the Reduced Representation genome
        -c pattern, --cut-site=pattern
                            Cutting sites of restriction enzyme. Ex: MspI(C-CGG),
                            Mael:(C-TAG), double-enzyme MspI&Mael:(C-CGG,C-TAG).
                            [Default: C-CGG]
        -L RRBS_LOW_BOUND, --low=RRBS_LOW_BOUND
                            Lower bound of fragment length (excluding C-CGG ends)
                            [Default: 20]
        -U RRBS_UP_BOUND, --up=RRBS_UP_BOUND
                            Upper bound of fragment length (excluding C-CGG ends)
                            [Default: 500]

      General options:
        -t TAG, --tag=TAG   [Y]es for undirectional lib, [N]o for directional
                            [Default: N]
        -s CUTNUMBER1, --start_base=CUTNUMBER1
                            The first cycle of the read to be mapped [Default: 1]
        -e CUTNUMBER2, --end_base=CUTNUMBER2
                            The last cycle of the read to be mapped [Default: 200]
        -a FILE, --adapter=FILE
                            Input text file of your adaptor sequences (to be
                            trimmed from the 3'end of the reads, ). Input one seq
                            for dir. lib., twon seqs for undir. lib. One line per
                            sequence. Only the first 10bp will be used
        --am=ADAPTER_MISMATCH
                            Number of mismatches allowed in adapter [Default: 0]
        -g GENOME, --genome=GENOME
                            Name of the reference genome (should be the same as
                            "-f" in bs_seeker2-build.py ) [ex. chr21_hg18.fa]
        -m NO_MISMATCHES, --mismatches=NO_MISMATCHES
                            Number of mismatches in one read [Default: 4]
        --aligner=ALIGNER   Aligner program for short reads mapping: bowtie,
                            bowtie2, soap, rmap [Default: bowtie]
        -p PATH, --path=PATH
                            Path to the aligner program. Detected:
                            bowtie: ~/install/bowtie
                            bowtie2: ~/install/bowtie2
                            rmap: ~/install/rmap/bin
                            soap: ~/install/soap/
        -d DBPATH, --db=DBPATH
                            Path to the reference genome library (generated in
                            preprocessing genome) [Default: ~/i
                            nstall/BSseeker2/bs_utils/reference_genomes]
        -l NO_SPLIT, --split_line=NO_SPLIT
                            Number of lines per split (the read file will be split
                            into small files for mapping. The result will be
                            merged. [Default: 4000000]
        -o OUTFILE, --output=OUTFILE
                            The name of output file [INFILE.bs(se|pe|rrbs)]
        -f FORMAT, --output-format=FORMAT
                            Output format: bam, sam, bs_seeker1 [Default: bam]
        --no-header         Suppress SAM header lines [Default: False]
        --temp_dir=PATH     The path to your temporary directory [Detected: /tmp]
        --XS=XS_FILTER      Filter definition for tag XS, format X,Y. X=0.8 and
                            y=5 indicate that for one read, if #(mCH sites)/#(all
                            CH sites)>0.8 and #(mCH sites)>5, then tag XS=1; or
                            else tag XS=0. [Default: 0.5,5]
        --multiple-hit      Output reads with multiple hits to
                            file"Multiple_hit.fa"
        -v, --version       show version of BS-Seeker2

      Aligner Options:
        You may specify any additional options for the aligner. You just have
        to prefix them with --bt- for bowtie, --bt2- for bowtie2, --soap- for
        soap, --rmap- for rmap, and BS Seeker will pass them on. For example:
        --bt-p 4 will increase the number of threads for bowtie to 4, --bt--
        tryhard will instruct bowtie to try as hard as possible to find valid
        alignments when they exist, and so on. Be sure that you know what you
        are doing when using these options! Also, we don't do any validation
        on the values.


##Examples :

* WGBS library ; alignment mode, bowtie ; map to WGBS index

        python bs_seeker2-align.py -i WGBS.fa --aligner=bowtie -o WGBS.bam -f bam -g genome.fa

* WGBS library ; alignment mode, bowtie2-local ; map to WGBS index

        python bs_seeker2-align.py -i WGBS.fa --aligner=bowtie2 -o WGBS.bam -f bam -g genome.fa

* WGBS library ; alignment mode, bowtie2-end-to-end ; map to WGBS index

        python bs_seeker2-align.py -i WGBS.fa -m 3 --aligner=bowtie2 -o WGBS.bam -f bam -g genome.fa --bt2--end-to-end

* RRBS library ; alignment mode, bowtie ; map to RR index

        python bs_seeker2-align.py -i RRBS.fa --aligner=bowtie -o RRBS.bam -g genome.fa -r -a adapter.txt

* RRBS library ; alignment mode, bowtie ; map to WG index

        python bs_seeker2-align.py -i RRBS.fa --aligner=bowtie -o RRBS.bam -g genome.fa -a adapter.txt

* RRBS library ; alignment mode, bowtie2-end-to-end ; map to WG index

        python bs_seeker2-align.py -i RRBS.fa --aligner=bowtie -o RRBS.bam -g genome.fa -a adapter.txt --bt2--end-to-end

* Align from qseq format for RRBS with bowtie, specifying lengths of fragments ranging [40bp, 400bp]

        python bs_seeker2-align.py -i RRBS.qseq --aligner=bowtie -o RRBS.bam -f bam -g genome.fa -r --low=40 --up=400 -a adapter.txt

The parameters '--low' and '--up' should be the same with corresponding parameters when building the genome index

* WGBS library ; alignment mode, bowtie ; map to WGBS index; use 8 threads for alignment

        python bs_seeker2-align.py -i WGBS.fa --aligner=bowtie -o WGBS.bam -f bam -g genome.fa --bt-p 4

BS-Seeker2 will run TWO bowtie instances in parallel.


##Input file:

- Adapter.txt (example) :

            AGATCGGAAGAGCACACGTC


##Output files:

- SAM file

    Sample:

        10918   0       chr1    133859922       255     100M    *       0       0       TGGTTGTTTTTGTTATAGTTTTTTGTTGTAGAGTTTTTTTTGGAAAGTTGTGTTTATTTTTTTTTTTGTTTGGGTTTTGTTTGAAAGGGGTGGATGAGTT        *       XO:Z:+FW        XS:i:0  NM:i:3  XM:Z:x--yx-zzzy--y--y--zz-zyx-yx-y--------z------------x--------z--zzz----y----y--x-zyx--------y--------z   XG:Z:-C_CGGCCGCCCCTGCTGCAGCCTCCCGCCGCAGAGTTTTCTTTGGAAAGTTGCGTTTATTTCTTCCCTTGTCTGGGCTGCGCCCGAAAGGGGCAGATGAGTC_AC


    Format descriptions:

        BS-Seeker2 specific tags:
        XO : orientation, from forward/reverted
        XS : 1 when read is recognized as not fully converted by bisulfite treatment, or else 0
        XM : number of sites for mismatch
                X: methylated CG
                x: un-methylated CG
                Y: methylated CHG
                y: un-methylated CHG
                Z: methylated CHH
                z: un-methylated CHH
        XG : genome sequences, with 2bp extended on both ends, from 5' to 3'
        YR : tag only for RRBS, serial id of mapped fragment
        YS : tag only for RRBS, start position of mapped fragment
        YE : tag only for RRBS, end position of mapped fragment

        Note:
            For reads mapped on Watson(minus) strand, the 10th colum in SAM file is not the original reads but the revered sequences.


##Tips:

- Removing adapter is recommended.

	If you don't know what's your parameter, please ask the person who generate the library for you.

	If you are too shy to ask for it, you can try to de novo motif finding tools (such as [DME](http://cb1.utdallas.edu/dme/index.htm) and [MEME](http://meme.nbcr.net/meme/cgi-bin/meme.cgi)) find the enriched pattern in 1000 reads.

	Of course, you can also use other tools (such as [cutadapt](http://code.google.com/p/cutadapt/) ) to remove adaptor first.

- It's always better to use a wider range for fragment length.

	For example, if 95% of reads come from fragments with length range [50bp, 250bp], you'd better choose [40bp, 300bp].

- Fewer mismatches for the 'local alignment' mode.

    As the 'local alignment', the bad sequenced bases are usually trimmed, and would not be considered by the parameter "-m".
    It is suggested to user fewer mismatches for the 'local alignment' mode.

(3) bs_seeker2-call_methylation.py
------------


This module calls methylation levels from the mapping result.

##Usage:

        $ python bs_seeker2-call_methylation.py -h
        Options:
          -h, --help            show this help message and exit
          -i INFILE, --input=INFILE
                                BAM output from bs_seeker2-align.py
          -d DBPATH, --db=DBPATH
                                Path to the reference genome library (generated in
                                preprocessing genome) [Default: ~/install
                                /BSseeker2/bs_utils/reference_genomes]
          -o OUTFILE, --output-prefix=OUTFILE
                                The output prefix to create ATCGmap and wiggle files
                                [INFILE]
          --wig=OUTFILE         The output .wig file [INFILE.wig]
          --CGmap=OUTFILE       The output .CGmap file [INFILE.CGmap]
          --ATCGmap=OUTFILE     The output .ATCGmap file [INFILE.ATCGmap]
          -x, --rm-SX           Removed reads with tag 'XS:i:1', which would be
                                considered as not fully converted by bisulfite
                                treatment [Default: False]
          --txt                 Show CGmap and ATCGmap in .gz [Default: False]
          -r READ_NO, --read-no=READ_NO
                                The least number of reads covering one site to be
                                shown in wig file [Default: 1]
          -v, --version         show version of BS-Seeker2


##Example :

-For WGBS (whole genome bisulfite sequencing):

        python bs_seeker2-call_methylation.py -i WGBS.bam -o output --db /path/to/BSseeker2/bs_utils/reference_genomes/genome.fa_bowtie/

-For RRBS:

        python bs_seeker2-call_methylation.py -i RRBS.bam -o output --db /path/to/BSseeker2/bs_utils/reference_genomes/genome.fa_rrbs_40_400_bowtie2/

-For RRBS and removed un-converted reads (with tag XS=1):

        python bs_seeker2-call_methylation.py -x -i RRBS.bam -o output --db /path/to/BSseeker2/bs_utils/reference_genomes/genome.fa_rrbs_75_280_bowtie2/

-For RRBS and only show sites covered by at least 10 reads in WIG file:

        python bs_seeker2-call_methylation.py -r 10 -i RRBS.bam -o output --db /path/to/BSseeker2/bs_utils/reference_genomes/genome.fa_rrbs_75_280_bowtie2/


The folder “genome.fa\_rrbs\_40\_500\_bowtie2” is built  in the first step

##Output files:

- wig file

        Sample:

        variableStep chrom=chr1
        3000419	0.000000
        3000423	-0.2
        3000440	0.000000
        3000588	0.5
        3000593	-0.000000


        Format descriptions:
        WIG file format. Negative value for 2nd column indicate a Cytosine on minus strand.


- CGmap file

    Sample:

        chr1	G	3000851	CHH	CC	0.1	1	10
        chr1	C	3001624	CHG	CA	0.0	0	9
        chr1	C	3001631	CG	CG	1.0	5	5

    Format descriptions:

        (1) chromosome
        (2) nucleotide on Watson (+) strand
        (3) position
        (4) context (CG/CHG/CHH)
        (5) dinucleotide-context (CA/CC/CG/CT)
        (6) methylation-level = #_of_C / (#_of_C + #_of_T).
        (7) #_of_C (methylated C, the count of reads showing C here)
        (8) = #_of_C + #_of_T (all Cytosines, the count of reads showing C or T here)


- ATCGmap file

    Sample:

        chr1	T	3009410	--	--	0	10	0	0	0	0	0	0	0	0	na
        chr1	C	3009411	CHH	CC	0	10	0	0	0	0	0	0	0	0	0.0
        chr1	C	3009412	CHG	CC	0	10	0	0	0	0	0	0	0	0	0.0
        chr1	C	3009413	CG	CG	0	10	50	0	0	0	0	0	0	0	0.83


    Format descriptions:

        (1) chromosome
        (2) nucleotide on Watson (+) strand
        (3) position
        (4) context (CG/CHG/CHH)
        (5) dinucleotide-context (CA/CC/CG/CT)

        (6) - (10) plus strand
        (6) # of reads from Watson strand mapped here, support A on Watson strand
        (7) # of reads from Watson strand mapped here, support T on Watson strand
        (8) # of reads from Watson strand mapped here, support C on Watson strand
        (9) # of reads from Watson strand mapped here, support G on Watson strand
        (10) # of reads from Watson strand mapped here, support N

        (11) - (15) minus strand
        (11) # of reads from Crick strand mapped here, support A on Watson strand and T on Crick strand
        (12) # of reads from Crick strand mapped here, support T on Watson strand and A on Crick strand
        (13) # of reads from Crick strand mapped here, support C on Watson strand and G on Crick strand
        (14) # of reads from Crick strand mapped here, support G on Watson strand and C on Crick strand
        (15) # of reads from Crick strand mapped here, support N

        (16) methylation_level = #C/(#C+#T) = C8/(C7+C8) for Watson strand, =C14/(C11+C14) for Crick strand;
        "nan" means none reads support C/T at this position.



Contact Information
============

If you still have questions on BS-Seeker 2, or you find bugs when using BS-Seeker 2, or you have suggestions, please write email to [Weilong Guo](guoweilong@gmail.com).



Questions & Answers
============

(1) Speed-up your alignment

Q: "It takes me days to do the alignment for one lane" ...

A: Yes, alignment is a time-consuming work, especially because the sequencing depth is increasing. An efficient way to align is :

    i. cut the original sequence file into multiple small pieces;

        Ex: split -l 4000000 input.fq

    ii. align them in parallel;
    iii. merge all the BAM files into a single one before running "bs-seeker2_call-methylation.py" (user "samtools merge" command).

        Ex: samtools merge out.bam in1.bam in2.bam in3.bam

(2) read in BAM/SAM

Q: Is the read sequence in BAM/SAM file is the same as my original one?

A: NO. They are different for several reasons.

    i. For RRBS, some reads are short because of trimming of the adapters
    ii. For read mapping on Crick (-) strand, the reads are in fact the antisense version of the original sequence, opposite both in nucleotides and direction

(3) "Pysam" package related problem

Q: I'm normal account user for Linux(Cluster). I can't install "pysam". I get following error massages:


        $ python setup.py install
        running install
        error: can't create or remove files in install directory
        The following error occurred while trying to add or remove files in the
        installation directory:
            [Errno 13] Permission denied: '/usr/lib64/python2.6/site-packages/test-easy-install-26802.write-test'
        ...


A: You can ask the administrator of your cluster to install pysam. If you don't want to bother him/her, you might need to build your own python, and then install the "pysam" package. The following script could be helpful for you.


        mkdir ~/install
        cd ~/install/

        # install python
        wget http://www.python.org/ftp/python/2.7.4/Python-2.7.4.tgz # download the python from websites
        tar zxvf Python-2.7.4.tgz # decompress
        cd Python-2.7.4
        ./configure --prefix=`pwd`
        make
        make install

        # Add the path of Python to $PATH
        #  Please add the following line to file ~/.bashrc

            export PATH=~/install/Python-2.7.4:$PATH

        # save the ~/.bashrc file
        source ~/.bashrc

        # install pysam package
        wget https://pysam.googlecode.com/files/pysam-0.7.4.tar.gz
        tar zxvf pysam-0.7.4.tar.gz
        cd pysam-0.7.4
        python setup.py build
        python setup.py install
        # re-login the shell after finish installing pysam

        # install BS-Seeker2
        wget https://github.com/BSSeeker/BSseeker2/archive/master.zip
        mv master BSSeeker2.zip
        unzip BSSeeker2.zip
        cd BSseeker2-master/


(4)Run BS-Seeker2

Q: Can I add the path of BS-Seeker2's *.py to the $PATH, so I can call
BS-Seeker2 from anywhere?

A: If you're using the "python" from path "/usr/bin/python", you can directly
add the path of BS-Seeker2 in file "~/.bash_profile" (bash) or "~/.profile"
(other shell) or "~/.bashrc" (per-interactive-shell startup).
But if you are using python under other directories, you might need to modify
BS-Seeker2's script first. For example, if your python path is "/my_python/python",
please change the first line in "bs_seeker-build.py", "bs_seeker-align.py" and
"bs_seeker-call_methylation.py" to

        #!/my_python/python

Then add

        export PATH=/path/to/BS-Seeker2/:$PATH

to file "~/.bash_profile" (e.g.), and source the file:

        source  ~/.bash_profile

Then you can use BS-Seeker2 globally by typing:

        bs_seeker_build.py -h
        bs_seeker-align.py -h
        bs_seeker-call_methylation.py -h

(5) Unique alignment

Q: If I want to only keep alignments that map uniquely, is this an argument I should supply directly
to Bowtie2 (via BS Seeker 2's command line option), or is this an option that's available in
BS Seeker 2 itself?

A: BS-Seeker2 reports unique alignment by default already. If you want to know how many reads
could have multiple hits, run BS-Seeker2 with parameter "--multiple-hit".


(6) Output

Q: In CGmap files, why some lines shown "--" but not a motif (CG/CHG/CHH), for example:

    chr01   C       4303711 --      CC      0.0     0       2
    chr01   C       4303712 --      CN      0.0     0       2

A: In this example, the site 4303713 is "N" in genome, thus we could not decide the explict pattern.

(7) Algorithm to remove the adapter.

Q: What's the algorithm to remove the adapter

A: BS-Seeker2 has built-in algorithm for removing the adapter, which is developed by [Weilong Guo](http://bioinfo.au.tsinghua.edu.cn/member/wguo/index.html).

    First, if the adapter was provided as "AGATCGGAAGAGCACACGTC", only the first 10bp would be used.
    Second, we use semi-local alignment strategy for removing the adapters.
    Exmaple:

        Read :       ACCGCGTTGATCGAGTACGTACGTGGGTC
        Adapter :    ....................ACGTGGGTCCCG

        no_mismatch : the maximum number allowed for mismatches

        Algorithm: (allowing 1 mismatch)
        -Step 1:
          ACCGCGTTGATCGAGTACGTACGTGGGTC
          ||XX
          ACGTGGGTCCCG
        -Step 2:
          ACCGCGTTGATCGAGTACGTACGTGGGTC
           X||X
          .ACGTGGGTCCCG
        -Step 3:
          ACCGCGTTGATCGAGTACGTACGTGGGTC
            XX
          ..ACGTGGGTCCCG
        -Step ...
        -Step N:
          ACCGCGTTGATCGAGTACGTACGTGGGTC
                              |||||||||
          ....................ACGTGGGTCCCG
        Success & return!

    Third, we also removed the synthesized bases at the end of RRBS fragments.
    Take the "C-CGG" cutting site as example,

        - - C|U G G - - =>cut=> - - C      =>add=> - - C|C G =>sequencing
        - - G G C|C - -         - - G G C          - - G G C

    In our algorithm, the "CG" in "--CCG" (upper strand) was trimmed, in order to get accurate methyaltion level.





