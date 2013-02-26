BS-Seeker 2
=========

BS-Seeker 2 performs accurate and fast mapping of bisulfite-treated short reads. BS-Seeker 2 is an updated version on BS-Seeker.

0. Availability
============

The source code for this package is available from
https://github.com/pfiziev/BS-Seeker. 
Also, you can use an instance of BS-Seeker 2 in Galaxy from http://galaxy.hoffman2.idre.ucla.edu. 
(Label: "NGS: Methylation Mapping"/"Methylation Map with BS Seeker2")


1. Remarkable new features
============
* Reduced index for RRBS, accelerating the mapping speed and increasing mappability
* Allowing local alignment with Bowtie 2, increased the mappability

2. Other features
============
* Supported library types
  - whole genomewide bisulfite sequencing (WGBS)
	- reduced representative bisulfite sequencing (RRBS)

* Supported formats for input file
	- [fasta](http://en.wikipedia.org/wiki/FASTA_format)
	- [fastq](http://en.wikipedia.org/wiki/FASTQ_format)
	- [qseq](http://jumpgate.caltech.edu/wiki/QSeq)
	- pure sequence

* Supported alignment tools
	- [bowtie](http://bowtie-bio.sourceforge.net/index.shtml)
	- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
	- [SOAP](http://soap.genomics.org.cn/)

* Supported fortmants for mapping results
	- [BAM](http://genome.ucsc.edu/FAQ/FAQformat.html#format5.1)
	- [SAM](http://samtools.sourceforge.net/)
	- [BS-seeker 1](http://pellegrini.mcdb.ucla.edu/BS_Seeker/USAGE.html)

3. System requirements
============

* Linux or Mac OS platform
* One of the following Aligner
  - [bowtie](http://bowtie-bio.sourceforge.net/)
  - [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/) (Recommend) 
  - [soap](http://soap.genomics.org.cn/)
* Python (Version 2.6 +)

  (It is normally pre-installed in Linux. Type " python -V" to see the installed version. Python is also freely available at http://www.python.org/download/)

* [pysam](http://code.google.com/p/pysam/) package needed. 
  


4. Modules' descriptions
============
(0) FilterReads.py 

Optional and independent module. 
Some reads would be extremely examplified during the PCR. This script helps you get unique reads before doing the mapping. You can decide whether or not to filter reads before doing the mapping.

	$ python FilterReads.py 
	Usage: FilterReads.py -i <input> -o <output> [-k]
	Author : Guo, Weilong; guoweilong@gmail.com; 2012-11-10
	Unique reads for qseq/fastq/fasta/sequencce, and filter 
	low quality file in qseq file.
	
	Options:
	  -h, --help  show this help message and exit
	  -i FILE     Name of the input qseq/fastq/fasta/sequence file
	  -o FILE     Name of the output file
	  -k          Would not filter low quality reads if specified


(1) bs_seeker2-build.py 

Module to build the index for BS-Seeker2. 
[Attention] Index built for BS-Seeker2 is different from the index for BS-Seeker 1.

For RRBS, you need to specify "-r" in the parameters. Also, you need to specify LOW_BOUND and UP_BOUND for the range of fragment lengths according your protocol. 
[Attention] The fragment length is different from read length. Fragments refers to the DNA fragements which you get by cutting from gel or other methods. Lengths of fragments are supposed to be in a range, such as [50bp,250bp].
[Attention] The indexes for RRBS and WGBS are different. Also, indexes for RRBS are variant for different values of parameters (LOW_BOUND and UP_BOUND).

	$ python bs_seeker2-build.py 
	
	Usage: bs_seeker2-build.py [options]
	
	Options:
	  -h, --help            show this help message and exit
	  -f FILE, --file=FILE  Input your reference genome file (fasta)
	  --aligner=ALIGNER     Aligner program to perform the analysis: bowtie,
	                        bowtie2, soap [bowtie2]
	  -p PATH, --path=PATH  Path to the aligner program. Defaults:
	                        bowtie: /home/guoweilong/bowtie-0.12.7/
	                        bowtie2: /home/guoweilong/bowtie-0.12.7/
	                        soap: /home/guoweilong/soap2.21release/
	  -d DBPATH, --db=DBPATH
	                        Path to the reference genome library (generated in
	                        preprocessing genome) [/home/guoweilong/BS-
	                        Seeker/bs_utils/reference_genomes]
	
	  Reduced Representation Bisulfite Sequencing Options:
	    Use this options with conjuction of -r [--rrbs]
	
	    -r, --rrbs          Preprocess the genome for analysis of Reduced
	                        Representation Bisulfite Sequencing experiments
	    -l LOW_BOUND, --low=LOW_BOUND
	                        lower bound [50]
	    -u UP_BOUND, --up=UP_BOUND
	                        upper bound [300]

Example

* Build Arabidoposis genome index for whole-genome bisulphite sequencing
        
        python bs_seeker2-build.py -f Arabidopsis.fa --aligner=bowtie2 -p ~/install/bowtie2-2.0.0-beta7/ 
        
* Build Arabidoposis genome index for RRBS library with default parameters
    
        python bs_seeker2-build.py -f Arabidopsis.fa --aligner=bowtie2 -p ~/install/bowtie2-2.0.0-beta7/ -r

* Build Arabidoposis genome index for RRBS library with fragment lengths ranging [50bp, 250bp]
    
        python bs_seeker2-build.py -f Arabidopsis.fa --aligner=bowtie2 -p ~/install/bowtie2-2.0.0-beta7/ -r -l 50 -u 250


(2) bs_seeker2-align.py 

	$ python bs_seeker2-align.py 
	
	Usage: bs_seeker2-align.py [options]
	
	Options:
	  -h, --help            show this help message and exit
	
	  For single end reads:
	    -i INFILE, --input=INFILE
	                        Input your read file name (FORMAT: sequences, illumina
	                        fastq, qseq,fasta)
	
	  For pair end reads:
	    -1 FILE, --input_1=FILE
	                        Input your read file end 1 (FORMAT: sequences,
	                        illumina fastq, qseq)
	    -2 FILE, --input_2=FILE
	                        Input your read file end 2 (FORMAT: sequences,
	                        illumina fastq, qseq)
	    --minins=MIN_INSERT_SIZE
	                        The minimum insert size for valid paired-end
	                        alignments [-1]
	    --maxins=MAX_INSERT_SIZE
	                        The maximum insert size for valid paired-end
	                        alignments [400]
	
	  Reduced Representation Bisulfite Sequencing Options:
	    -r, --rrbs          Process reads from Reduced Representation Bisulfite
	                        Sequencing experiments
	    --rrbs-tag=TAG      Msp-I tag: CGG TGG CGA or CGG/TGG (both)
	    --low=RRBS_LOW_BOUND
	                        lower bound [50]
	    --up=RRBS_UP_BOUND  upper bound [300]
	
	  General options:
	    -t TAG, --tag=TAG   Yes for undirectional lib, no for directional [Y]
	    -s CUTNUMBER1, --start_base=CUTNUMBER1
	                        The first base of your read to be mapped [1]
	    -e CUTNUMBER2, --end_base=CUTNUMBER2
	                        The last cycle number of your read to be mapped [200]
	    -a FILE, --adapter=FILE
	                        Input text file of your adaptor sequences (to be
	                        trimed from the 3'end of the reads). Input 1 seq for
	                        dir. lib., 2 seqs for undir. lib. One line per
	                        sequence
	    -g GENOME, --genome=GENOME
	                        Name of the reference genome (the same as the
	                        reference genome file in the preprocessing step) [ex.
	                        chr21_hg18.fa]
	    -m INT_NO_MISMATCHES, --mismatches=INT_NO_MISMATCHES
	                        Number of mismatches (0,1,...,read length) [4]
	    --aligner=ALIGNER   Aligner program to perform the analisys: bowtie,
	                        bowtie2, soap [bowtie2]
	    -p PATH, --path=PATH
	                        Path to the aligner program. Defaults:
	                        bowtie: /home/guoweilong/bowtie-0.12.7/
	                        bowtie2: /home/guoweilong/bowtie-0.12.7/
	                        soap: /home/guoweilong/soap2.21release/
	    -d DBPATH, --db=DBPATH
	                        Path to the reference genome library (generated in
	                        preprocessing genome) [/home/guoweilong/BS-
	                        Seeker/bs_utils/reference_genomes]
	    -l NO_SPLIT, --split_line=NO_SPLIT
	                        Number of lines per split (the read file will be split
	                        into small files for mapping. The result will be
	                        merged. [4000000]
	    -o OUTFILE, --output=OUTFILE
	                        The name of output file [INFILE.bs(se|pe|rrbs)]
	    -f FORMAT, --output-format=FORMAT
	                        Output format: bam, sam, bs_seeker1 [bam]
	    --no-header         Suppress SAM header lines [False]
	    --temp_dir=PATH     The path to your temporary directory [/tmp]
	    --XS=XS_FILTER      Filter definition for tag XS, format X,Y. X=0.8 and
		                y=5 indicate that for one read, if #(mCH sites)/#(all
		                CH sites)>0.8 and #(mCH sites)>5, then tag XS=1; or
		                else tag XS=0. [0.5,5]
	
	  Aligner Options:
	    You may specify any additional options for the aligner. You just have
	    to prefix them with --bt- for bowtie, --bt2- for bowtie2, --soap- for
	    soap, and BS Seeker will pass them on. For example: --bt-p 4 will
	    increase the number of threads for bowtie to 4, --bt--tryhard will
	    instruct bowtie to try as hard as possible to find valid alignments
	    when they exist, and so on. Be sure that you know what you are doing
	    when using these options! Also, we don't do any validation on the
	    values.

Examples

Align from fasta format with bowtie2 for whole genome, allowing 3 mismatches

        python bs_seeker2-align.py -i input.fa -m 3 --aligner=bowtie2 -p ~/install/bowtie2/ -o output.bam -f bam -g Arabidopsis.fa 

Align from qseq format for RRBS, allowing 5 mismatches, default parameters for RRBS fragment

        python bs_seeker2-align.py -i input.qseq -m 5 --aligner=bowtie2 -p ~/install/bowtie2/ -o output.bam -f bam -g Arabidopsis.fa -r

Align from qseq format for RRBS, allowing 5 mismatches, specifying lengths of frament

        python bs_seeker2-align.py -i input.qseq -m 5 --aligner=bowtie2 -p ~/install/bowtie2/ -o output.bam -f bam -g Arabidopsis.fa -r --low=50 --up=250

The parameters '--low' and '--up' should be the same with correponding parameters when building the genome index


- SAM file
Sample:
	
        10918   0       chr1    133859922       255     100M    *       0       0       TGGTTGTTTTTGTTATAGTTTTTTGTTGTAGAGTTTTTTTTGGAAAGTTGTGTTTATTTTTTTTTTTGTTTGGGTTTTGTTTGAAAGGGGTGGATGAGTT        *       XO:Z:+FW        XS:i:0  NM:i:3  XM:Z:x--yx-zzzy--y--y--zz-zyx-yx-y--------z------------x--------z--zzz----y----y--x-zyx--------y--------z   XG:Z:-C_CGGCCGCCCCTGCTGCAGCCTCCCGCCGCAGAGTTTTCTTTGGAAAGTTGCGTTTATTTCTTCCCTTGTCTGGGCTGCGCCCGAAAGGGGCAGATGAGTC_AC


Format:

        BSseeker2 specific tags:
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

(3) bs_seeker2-call_methylation.py

This module calls methylation levels from the mapping result.

	$ python bs_seeker2-call_methylation.py 
	
	Usage: bs_seeker2-call_methylation.py [options]
	
	Options:
	  -h, --help            show this help message and exit
	  -i INFILE, --input=INFILE
	                        BAM output from bs_seeker2.py
	  --db=DBPATH           Path to the reference genome library (generated in
	                        preprocessing genome) [/home/guoweilong/BS-
	                        Seeker/bs_utils/reference_genomes]
	  -o OUTFILE, --output-prefix=OUTFILE
	                        The output prefix to create ATCGmap and wiggle files
	                        [INFILE]
	  --wig=OUTFILE         The output .wig file [INFILE.wig]
	  --CGmap=OUTFILE       The output .CGmap file [INFILE.CGmap]
	  --ATCGmap=OUTFILE     The output .ATCGmap file [INFILE.ATCGmap]
	  -x, --rm-SX           Removed reads with tag 'XS:i:1', which would be
		                considered as not fully converted by bisulfite
		                treatment
	  -r READ_NO, --read-no=READ_NO
		                The least number of reads covering one site to be
		                shown in wig file [1]

Example :

For WGBS (whole genome bisulfite sequencing):

        python bs_seeker2-call_methylation.py -i ath_whole.bam -o Ath_whole --db /path/to/BSseeker2/bs_utils/reference_genomes/Arabidopsis.fa_bowtie2/

For RRBS:

        python bs_seeker2-call_methylation.py -i ath_rrbs.bam -o output --db /path/to/BSseeker2/bs_utils/reference_genomes/Arabidopsis.fa_rrbs_75_280_bowtie2/ 

For RRBS and removed un-converted reads (with tag XS=1):

        python bs_seeker2-call_methylation.py -x -i ath_rrbs.bam -o output --db /path/to/BSseeker2/bs_utils/reference_genomes/Arabidopsis.fa_rrbs_75_280_bowtie2/ 

For RRBS and only show sites covered by at least 10 reads in WIG file:

        python bs_seeker2-call_methylation.py -r 10 -i ath_rrbs.bam -o output --db /path/to/BSseeker2/bs_utils/reference_genomes/Arabidopsis.fa_rrbs_75_280_bowtie2/ 


The folder “Arabidopsis.fa_rrbs_75_280_bowtie2” is builded  in the first step

Description of output files:

- wig file

        variableStep chrom=chr1
        3000419	0.000000
        3000423	-0.2
        3000440	0.000000
        3000588	0.5
        3000593	-0.000000


Format:
WIG file format. Negative value for 2nd column indicate a Cytosine on minus strand.


- CGmap file
Sample:
	
        chr1	G	3000851	CHH	CC	0.1	1	10
        chr1	C	3001624	CHG	CA	0.0	0	9
        chr1	C	3001631	CG	CG	1.0	5	5
	
Format description:

        (1) chromosome 
        (2) nucleotide 
        (3) position 
        (4) context 
        (5) dinucleotide-context 
        (6) methyltion-level = #-of-C / (#-of-C + #-of-T)
        (7) #-of-C (methylated)
        (8) (#-ofC + #-of-T) (all cytosines)


- ATCGmap file
Sample:
        
        chr1	T	3009410	--	--	0	1	0	0	0	0	0	0	0	0	na
        chr1	C	3009411	CHH	CC	0	1	0	0	0	0	0	0	0	0	0.0
        chr1	C	3009412	CHG	CC	0	1	0	0	0	0	0	0	0	0	0.0
        chr1	C	3009413	CG	CG	0	1	5	0	0	0	0	0	0	0	0.833333333333


Format description:

        (1) chromosome
        (2) nucleotide
        (3) position
        (4) context
        (5) dinucleotide-context

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

        (16) methylation_level = #C/(#C+#T); "nan" means none reads support C/T at this position.



Contact Information
============

If you still have questions on BS-Seeker 2, please write email to guoweilong@gmail.com.



