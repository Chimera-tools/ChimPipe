.. _manual:

======
Manual
======

This section describes the files ChimPipe takes as input, how to generate them, run the pipeline and interpret its output. 

Input files
===========
ChimPipe needs 4 mandatory types of files plus 1 optional.  

Mandatory:

* Paired-end (PE) RNA-seq reads
* Genome index 
* Genome annotation
* Transcriptome annotation index

Optional:

* Similarity between gene pairs text file


Paired-end (PE) RNA-seq reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ChimPipe has been designed to deal with `Illumina paired-end`_ RNA sequencing data. It takes two `FASTQ`_ files as input, one for the first mates and another one for the second mates in the read pairs respectively. Make sure that they are in the same directory and they are named according to this convention: 

.. _Illumina paired-end: http://technology.illumina.com/technology/next-generation-sequencing/paired-end-sequencing_assay.ilmn
.. _FASTQ: http://maq.sourceforge.net/fastq.shtml

* **Mate 1**. "SampleId + [.1|_1] + [.fastq|.txt] (+ .gz if compressed)"
* **Mate 2**. "SampleId + [.2|_2] + [.fastq|.txt] (+ .gz if compressed)"

E. g. BERGER_1.fastq.gz and BERGER_2.fastq.gz would be the compressed FASTQ files for mate 1 and mate 2 in the BERGER sample. Note that you should use the same convention for both mates, so if mate 1 is BERGER_1.fastq.gz, mate 2 can not be BERGER.2.txt.gz

The FASTQ file uses four lines per sequencing read. You need to check the format of the first line of each read, which begins with the '@' character and is followed by a read identifier. This identifier should meet one of the two Illumina standards to specify which member of the pair the read is:

* Illumina CASAVA package lower than 1.8. The identifier has to be a single string ended in /1 or /2 for mate 1 and mate 2, respectively. E. g.:

.. code-block:: bash
	
	$ # Mate 1
	@SRR018259.1_BI:080831_SL-XAN_0004_30BV1AAXX:5:1:708:1868#0/1
	GTAACATATTCACAGACATGTCAGTGTGCACTCAGGAACACTGGTTTCATT
	+SRR018259.1_BI:080831_SL-XAN_0004_30BV1AAXX:5:1:708:1868#0/1
	IIIIIIIIIIIIIII>=CII=8=H032/-++D+'.@)2+4/+)1'4.#"*.
	
	$ # Mate 2
	@SRR018259.1_BI:080831_SL-XAN_0004_30BV1AAXX:5:1:708:1868#0/2
	CAGATGGCATTGCAGAACTAAAAAAGCAAGCTGAATCAGTGTTAGATCTCC
	+SRR018259.1_BI:080831_SL-XAN_0004_30BV1AAXX:5:1:708:1868#0/2
	IIIGIIIIIIIIIFI:>DID<AFH=>78@I41I055549.?+.42-'%**'
	

* Illumina CASAVA package 1.8 and higher. The identifier consists on two strings separated by a space with the second one specifying the mate in the first digit. E. g:   

.. code-block:: bash
	
	# Mate 1
	@seq.1_set1:140:D034KACXX:1:2105:18345:62400 1:N:0:
	CCCAGCCTGTTTCCCTGCCTGGGAAACTAGAAAGAGGGCTTCTTTCTCCT
	+
	IJJJIIIIIGIGIEBHHHGFFFF=CDEEEEEDDDDCCCDDA?55><CBBB
	
	# Mate 2
	@seq.1_set1:140:D034KACXX:1:2105:18345:62400 2:N:0:
	GCACCCTTCACTCCCTCCCTTGGGCGCCTCCCTCCCGAGGGTAGGGACCC
	+
	FFHIJJCHIIHBHIIIAHFFFFFCDEDEECDBB;??@CD?CCCCCCC@CC

I case your FASTQ files do not meet this format you should modify the identifier. Awk is a perfect tool for such kind of problems. E. g: I downloaded a dataset from the NCBI Sequence Read archive:

.. code-block:: bash
	

	# Mate 1 read without a proper identifier. It has three strings as identifier and does not end with "/1"
	
	@SRR018259.1 BI:080831_SL-XAN_0004_30BV1AAXX:5:1:708:1868 length=51
	GTAACATATTCACAGACATGTCAGTGTGCACTCAGGAACACTGGTTTCATT
	+SRR018259.1 BI:080831_SL-XAN_0004_30BV1AAXX:5:1:708:1868 length=51
	IIIIIIIIIIIIIII>=CII=8=H032/-++D+'.@)2+4/+)1'4.#"*.
	
	# No worries, I can use awk to fix it. 
	
	$ awk '{if (NR%2==1){print $1"_"$2"_"$3"#0/1"} else {print $0}} dataset_1.fastq 		
	
	@SRR018259.1_BI:080831_SL-XAN_0004_30BV1AAXX:5:1:708:1868_length=51#0/1
	GTAACATATTCACAGACATGTCAGTGTGCACTCAGGAACACTGGTTTCATT
	+SRR018259.1_BI:080831_SL-XAN_0004_30BV1AAXX:5:1:708:1868_length=51#0/1
	IIIIIIIIIIIIIII>=CII=8=H032/-++D+'.@)2+4/+)1'4.#"*.

	$ # Finally, I apply the same procedure for the mate 2..

Genome index
~~~~~~~~~~~~
An indexed reference genome in GEM format has to be provided to do the mapping steps. You just need to run the *GEMtools indexer* (supplied with ChimPipe) with your genome in FASTA format to produce it:

.. _FASTA:
 
.. code-block:: bash

	$ gemtools index -i genome.fa 

It will produce 3 files in the directory where the genome is placed:

* **genome.gem** – indexed genome in GEM format (needed for running ChimPipe).   
* genome.hash – hash table with the genome (no needed). 
* genome.log – indexer log file.    

**Note:** It is recommended to use multiple threads with the option -t. 


Genome annotation
~~~~~~~~~~~~~~~~~~
Chimpipe also takes as input a genome annotation in `GTF`_ format with the annotated exons. It can contain other features different from exons, i. e. introns or UTR, but they will be not considered by the pipeline in the chimera detection process. This annotation has to contain at least one tag-value pair in the attributes field with the gene id and two optional pairs will be taken into account by ChimPipe if supplied: gene name and gene type. E.g:

.. _GTF: http://www.ensembl.org/info/website/upload/gff.html

.. code-block:: bash
	
	# This is an example of one annotated exon with an appropiated format. 	
	# The attributes are the gene id (mandatory), the gene type and gene name (optional), 
	# plus some additional tag-value pairs that will not be considered by ChimPipe.   
	
	chr1	HAVANA	exon	69091	70008	.	+	.	gene_id "ENSG00000186092.4"; transcript_id "ENST00000335137.3"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "OR4F5";
	transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "OR4F5-001"; exon_number 1; exon_id "ENSE00002319515.1"; level 2; tag "basic"; tag "appris_principal";
	tag	"CCDS"; ccdsid "CCDS30547.1"; havana_gene "OTTHUMG00000001094.1"; havana_transcript "OTTHUMT00000003223.1";

ChimPipe has been benchmarked with `Gencode v10`_ and `UCSC Known Genes`_ annotation. It displayed a better sensitivity with Gencode v10 while the similar false positive rate was similar (see Benchmark section). Thus, we encourage the user to use Gencode annotation, it is a richer annotation what increase the sensitivity of the chimera detection process. 

.. _Gencode v10: http://www.gencodegenes.org/releases/10.html
.. _UCSC Known Genes: https://genome.ucsc.edu/cgi-bin/hgTables?command=start

Transcriptome annotation index
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
An indexed transcriptome annotation in GEM format has to be given as input to find reads spanning annotated splice junctions. You only have to run the *GEMtools transcriptome indexer* (supplied with ChimPipe) with your previously generated GEM windexed genome and its annotation in GTF format to generate it. 

.. code-block:: bash

	$ $gemtools t-index -i genome.gem -a annotation.gtf	

It will produce 5 files in your current working directory:

* annotation.gtf.junctions – annotated splice junctions coordinates (no needed)
* annotation.gtf.junctions.fa – annotated splice junctions sequence (no needed)
* **annotation.gtf.junctions.gem** – transcriptome index in GEM format (needed)
* **annotation.gtf.junctions.keys** – keys to convert from transcriptome to genome (needed)
* annotation.gtf.junctions.log – indexer log file

**Note**: It is recommended to use multiple threads with the option -t. 

**IMPORTANT**: The indexed transcriptome annotation has to be placed in the same folder as the genome annotation to be used by ChimPipe.

Similarity between gene pairs (Optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Execute pipeline
================

1. Set up the environment
~~~~~~~~~~~~~~~~~~~~~~~~~
As explained in the :ref:`installation` section, you need to have installed BEDtools and SAMtools to execute ChimPipe, plus blast in case you want to produce your own similarity between transcript pairs text files (See **Similarity between gene pairs**). In case you do not have them, you can download an install them from their webpages. Once installed, you have to export the path to their binaries. 

Please check our :ref:`FAQ` section in case you have any problem to do it.  

2. Check the quality offset in your dataset   
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The quality scores (Q) measure the probability that a base is called incorrectly by the sequencing machine. Within your FASTQ files, they are represented in the fourth line of each read as an string of ASCII characters (each character correspond to the Q score of a certain base in the sequencing read). The correspondence between each ASCII character and the Q score is based on some offset. These offset vary with the sequencing platform (current Illumina machines uses 33, while older ones 33). 

**NOTE**: ChimPipe needs to know the offset used in your RNA-seq dataset to do the mapping steps. If you do not have this information, we provide a short script to easily check it (see :ref:`FAQ` section). 

3. Check the RNA-seq library type
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Different protocols that can be used to generate a RNA-seq library. There are also important differences among them that have to be taken into account in several steps of the chimera detection pipeline. However, ChimPipe can not determine the protocol used to produce your reads, so you need to supply this information manually with the option **--read-directionality <STRING>**. Where **STRING** has to be one of these library types:

* **NONE**. Not strand-specific protocol (unstranded data). The information about from which strand the transcript is transcribed is not available. Default configuration.

Strand-specific protocols (stranded data):
 
* **SENSE**. Transcript directly sequenced. Reads map to the transcript strand.
* **ANTISENSE**. Reverse complementary sequence of the transcript sequenced. Reads map to the opposite strand of the transcript. 
* **MATE1_SENSE**. Reads on the left of the fragment (mates 1) sequenced from the transcript (map in the transcrip strand), and the ones in the right (mates 2) sequenced from the complementary reverse sequence (map in the opposite strand). 
* **MATE2_SENSE**. Reads on the left of the fragment (mates 1) sequenced from the complementary reverse sequence (map in the opposite strand), and the ones in the right (mates 2) sequenced from the transcript (map in the transcrip strand). 
	
**NOTE**: if you do not know, you can ask your RNA-seq data provider or use our bash script (see :ref:`FAQ` section).
	
4. Run ChimPipe
~~~~~~~~~~~~~~~

.. code-block:: bash

	** Mandatory arguments:
		-i|--input			<INPUT_FILE>	First mate sequencing reads in FASTQ format. Pipeline designed to deal with paired-end
	 						data. Please make sure the second mate file is in the same directory as the first mate 
	 						file, and the files are named "YourSampleId_1.fastq.gz" and "YourSampleId_2.fastq.gz" 
	 						respectively. YourSampleId has to be provided with the -e argument.
		-g|--genome-index		<GEM>		Index for the reference genome in GEM format.
		-a|--annotation			<GTF>		Reference gene annotation file in GTF format.
		-q|--quality			<NUMBER>	Quality offset of the FASTQ files [33 | 64 | ignore].
		-e|--sample-id			<STRING>	Sample identifier (the output files will be named according to this id).    
		
	** [OPTIONS] can be:
	Reads information:
		-s|--stranded			<FLAG>		Flag to specify whether data is stranded. Default false (unstranded).
		--read-directionality		<STRING>	Directionality of the reads [MATE1_SENSE | MATE2_SENSE | MATE_STRAND_CSHL | SENSE | ANTISENSE | NONE]. Default NONE.
		--max-read-length		<NUMBER>	Maximum read length. This is used to create the de-novo transcriptome and acts as an upper bound. Default 150.
		
	Mapping parameters:
		-M|--mism-contiguous-map	<NUMBER>	Maximum number of mismatches for the contiguous mapping steps with the GEM mapper. Default 4?. Not working
		-m|--mism-split-map		<NUMBER>	Maximum number of mismatches for the segmental mapping steps with the GEM rna-mapper. Default 4?.	Not working
		-c|--consensus-splice-sites	<(couple_1)>, ... ,<(couple_s)>	with <couple> := <donor_consensus>+<acceptor_consensus>
	                                 			(list of couples of donor/acceptor splice site consensus sequences, default='(GT,AG),(GC,AG),(ATATC,A.),(GTATC,AT)'
		--min-split-size		<NUMBER>	Minimum split size for the segmental mapping steps. Default 15.
		--stats				<FLAG>		Enable mapping statistics. Default disabled.
		
	Chimeric junctions filter:
			--filter-chimeras		<STRING>	Configuration for the filtering module. Quoted string with 4 numbers separated by commas and ended in semicolom, 
								i.e. "1,2,75:50;", where:
												
									1st: minimum number of staggered reads spanning the chimeric junction.
									2nd: minimum number of paired-end reads encompassing the chimeric junction.		
									3rd: maximum similarity between the connected genes.
									4rd: maximum length of the high similar region between the connected genes.
		
								All these conditions have to be fulfilled for a chimeric junction to pass the filter. It is also possible to make 
								complex condifions by setting two different conditions where at least one of them has to be fulfilled. 
								I.e "10,0,0:0;1,1,0:0;". Default "5,0,80:30;1,1,80:30;".	
			--similarity-gene-pairs	<TEXT>			Text file containing similarity information between the gene pairs in the annotation. Needed for the filtering module 
								to discard junctions connecting highly similar genes. If not provided the junctions will not be filtered according to this criteria. 

	General:
		-o|--output-dir			<PATH>		Output directory. Default current working directory.
		--tmp-dir			<PATH>		Temporary directory. Default /tmp.
		-t|--threads			<PATH>		Number of threads to use. Default 1.
		-l|--log			<PATH>		Log level [error |warn | info | debug). Default info.
		--dry				<FLAG>		Test the pipeline. Writes the command to the standard output.
		--help				<FLAG>		Display usage information.




Output
======

By default, ChimPipe produces 3 files as output:

* First mapping BAM file
* Second mapping MAP file
* Chimeric junctions text file

**NOTE**: you can use the optional flag option **--no-cleanup** to ouput also intermediate files. 

First mapping BAM file
~~~~~~~~~~~~~~~~~~~~~~
`BAM`_ file containing the reads mapped in the genome, transcriptome and *de novo* transcriptome with the *GEMtools RNA pipeline*. 

Many next-generation sequencing analysis tools work with this format, so it can be used to do very different analyses such as gene and transcript quantification or differential gene expression analysis.

.. _BAM: http://samtools.github.io/hts-specs/SAMv1.pdf

Second mapping MAP file
~~~~~~~~~~~~~~~~~~~~~~~
MAP file containing reads segmentally mapped in the genome allowing for interchromosomal, different strand and unexpected genomic order mappings. 


Chimeric junctions text file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Tabular text file containing the detected chimeric junctions in your RNA-seq dataset. It has rows of 19 fields, where each row corresponds to a chimeric junction and the fields contains information about it. Here is a brief description of the 19 fields:

1. **juncId** - Chimeric junction identifier. It is an string encoding the position of the chimeric junction in the genome as follows: chrA"_"breakpointA"_"strandA":"chrB"_"breakpointB"_"strandB. E. g., "chr4_90653092_+:chr17_22023757_+" is a chimeric junction between the position 90653092 of the chromosome 4 in the plus strand, and the position 22023757 of the chromosome chr17 in the plus strand. 
2. **nbstag** - Number of staggered reads supporting the chimera.
3. **nbtotal** - Total number of reads supporting the chimera.
4. **maxbeg** - Maximum beginning of the chimeric junction,  The starting position at which 
5. **maxEnd** - Maximum end of the junction
6. **samechr** - Flag to specify if the connected gene pairs are in the same cromosome (1) or not (0).
7. **samestr** - Flag to specify if the connected gene pairs are in the same strand (1) or not (0), NA in case the *samechr* field was 0.
8. **dist** - Distance between the two breakpoints, NA in case the "samestr" field was 0.
9. **ss1** - Splice donor site sequence.
10. **ss2**	- Splice acceptor site sequence.
11. **gnlist1** - List of genes overlapping the first part of the chimera. 	
12. **gnlist2**	- List of genes overlapping the second part of the chimera. 
13. **gnname1** - Name of the genes in the field *gnlist1*, "." if unknown. 
14. **gnname2**	- Name of the genes in the field *gnlist1*, "." if unknown.
15. **bt1** - Biotype of the genes in the field *gnlist1*, "." if unknown. 
16. **bt2**	- Biotype of the genes in the field *gnlist2*, "." if unknown.
17. **PEsupport** - Total number of read pairs supporting the chimera, "." if not Paired-end support. It is a string containing information about the number of read pairs supporting the connection between the involved gene pairs as follows: geneA1-GeneA2:nbReadPairs,geneB1-geneB2:nbReadPairs. E.g.: "1-1:1,3-1:2" means that the connection between the genes 1, in the *gnlist1* and *gnlist2* respectively, is supported by 1 read pair; and the connection between the gene 3 in the *gnlist1* and the gene 1 in the *gnlist2* is supported by 2 read pairs. 
18. **maxSim** - Maximum percent of similarity in the BLAST alignment between the transcript with the longest BLAST alignment, "." if no blast hit found.
19. **maxLgal** - Maximum length of the BLAST alignment between all the transcripts of the gene pairs connected by the chimeric junction, "." if no blast hit found. 

**Example**

.. code-block:: bash

	Here is an example of a chimeric junction detected by ChimPipe

	juncId	nbstag	nbtotal	maxbeg	maxEnd	samechr	samestr	dist	ss1	ss2	gnlist1	gnlist2	gnname1	gnname2	bt1	bt2	PEsupport	maxSim	maxLgal
	chr1_121115975_+:chr1_206566046_+ 1 1 121115953 206566073 1 1 85450071 GC AG SRGAP2D, SRGAP2,SRGAP2C, SRGAP2D, SRGAP2,SRGAP2C, . . 1-1:2,1-2:2, 99.44 1067

