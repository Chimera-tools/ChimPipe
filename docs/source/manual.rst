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
ChimPipe has been designed to deal with `Illumina paired-end`_ RNA sequencing data. It takes two `FASTQ`_ files as input, one for the first mates and another one for the second mates in the read pairs respectively. They have to be located in the same directory and named according to this convention: 

.. _Illumina paired-end: http://technology.illumina.com/technology/next-generation-sequencing/paired-end-sequencing_assay.ilmn
.. _FASTQ: http://maq.sourceforge.net/fastq.shtml

* **Mate 1**. "SampleId + [.1|_1] + [.fastq|.txt] (+ .gz if compressed)"
* **Mate 2**. "SampleId + [.2|_2] + [.fastq|.txt] (+ .gz if compressed)"

E. g. BERGER_1.fastq.gz and BERGER_2.fastq.gz would be the compressed FASTQ files for mate 1 and mate 2 in the BERGER sample. 

.. warning:: Please make sure the 2 FASTQ files are in the **same directory** and you use the **same convention** for both mates. E. g: if mate 1 is BERGER_1.fastq.gz, mate 2 can not be BERGER.2.txt.gz

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

We also provide some **pre-generated genome indices** for human, mouse and drosophila genomes in the :ref:`Downloads` section. 

.. tip:: If your machine has more than one CPU it is recommended to run the indexer with multiple threads. Use the option ``-t <threads>``, where **threads** is the number of CPUs available. 


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

.. note:: ChimPipe has been benchmarked with `Gencode v10`_ and `UCSC Known Genes`_  humna annotation. It displayed a better sensitivity with Gencode v10 while the similar false positive rate was similar. Thus, it is advisable to use Gencode annotation, it is a richer annotation what increase the sensitivity of the chimera detection process. 

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

We also provide some **pre-generated transcriptome indices** for human, mouse and drosophila annotations in the :ref:`Downloads` section. 

.. tip:: If your machine has more than one CPU it is recommended to run the indexer with multiple threads. Use the option ``-t <threads>``, where **threads** is the number of CPUs available. 

.. warning:: The indexed transcriptome annotation has to be placed in the same folder as the genome annotation to be used by ChimPipe.

Similarity between gene pairs (Optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
One of ChimPipe's filtering steps to discard actefactual chimeras is to filter out those chimeric junctions connecting genes that encode transcripts with a high sequence homology. Although it is an optional filter, it is **strongly recommended**, since accoding our benchmark it improves much the specificity with a minimal decrease of the sensitivity.  

To enable this homology-based filtering you only need to run ChimPipe with the option ``--similarity-gene-pairs <TEXT FILE>``, where **TEXT FILE** is a file containing the matrix with information about the sequence similarity between gene pairs. 

You can download our **pre-generated matrices** por human, mouse and drosophila annotations from :ref:`Downloads` section or you can produce your own matrix with the script ``ChimPipe/src/bash/tools/similarity_bt_gnpairs.sh`` as follows:

	$ bash similarity_bt_gnpairs.sh annot.gtf genome.gem

This script will produce the matrix through 4 steps:

1. Extract the cDNA sequence of each transcript in the annotation.

2. Make a BLAST database out of the transcript sequences. 

3. Run BLAST on all trancript against all transcripts to detect local similarity between transcripts.

4. Produce a 8 fields matrix where each row corresponds to a gene pair and it contains information about the alignment between the pair of transcripts of this two genes with the maximum alignment similarity and length. Here is a brief description of the 8 fields:

	1. Gene name A
	2. Gene name B
	3. Transcripts alignment similarity
	4. Transcript alignment length
	5. Transcript name A
	6. Transcript name B
	7. Trancript A exonic length
	8. Transcript B exonic length

**Example** 

ENSG00000000003.10 ENSG00000003402.15 91.43 70 ENST00000373020.4 ENST00000309955.3 2206 14672
	
.. warning:: Make sure you run ChimPipe with a similarity matrix generated from the annotation and genome you are using.  

Execute pipeline
================

1. Set up the environment
~~~~~~~~~~~~~~~~~~~~~~~~~
As explained in the :ref:`installation` section, you need to have installed BEDtools and SAMtools to execute ChimPipe, plus blast in case you want to produce your own similarity between gene pairs text files (See **Similarity between gene pairs**). In case you do not have them, you can download an install them from their webpages. Once installed, you have to export the path to their binaries. 

Please check our :ref:`FAQ` section in case you have any problem.  

2. Check the quality offset in your dataset   
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The quality scores (Q) measure the probability that a base is called incorrectly by the sequencing machine. Within your FASTQ files, they are represented in the fourth line of each read as an string of ASCII characters (each character correspond to the Q score of a certain base in the sequencing read). The correspondence between each ASCII character and the Q score is based on some offset. These offset vary with the sequencing platform (current Illumina machines uses 33, while older 33). 

.. tip:: ChimPipe needs to know the offset considered in your RNA-seq data to do the mapping steps. If you do not have this information, a short script is provided to easily check it (see :ref:`FAQ` section). 

3. Check the RNA-seq library type
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Different protocols that can be used to generate a RNA-seq library. There are also important differences among them that have to be taken into account in several steps of the chimera detection pipeline. However, ChimPipe can not determine the protocol used to produce your reads, so you need to supply this information manually with the option ``--read-directionality <library>``. Where **library** has to be one of these library types:

* **NONE**. Not strand-specific protocol (unstranded data). The information about from which strand the transcript is transcribed is not available. Default configuration.

Strand-specific protocols (stranded data):
 
* **MATE1_SENSE**. Reads on the left of the fragment (mates 1) sequenced from the transcript (map in the transcrip strand), and the ones in the right (mates 2) sequenced from the complementary reverse sequence (map in the opposite strand). 
* **MATE2_SENSE**. Reads on the left of the fragment (mates 1) sequenced from the complementary reverse sequence (map in the opposite strand), and the ones in the right (mates 2) sequenced from the transcript (map in the transcrip strand). 
	
.. tip:: In case you do not know the type of library, use the bash script provided along ChimPipe (see :ref:`FAQ` section) or ask your RNA-seq data provider.
	
4. Run ChimPipe
~~~~~~~~~~~~~~~
Once you have the genome and transcriptome indices prepared, and you know the quality offset and the library type of your PE RNA-seq reads you can run ChimPipe as follows:

.. code-block:: bash
	
	bash ChimPipe.sh -i reads_1.fastq -g genome.gem -a annotation.gtf -q 33 -l UNSTRANDED -e sample1

All these files and parameters given as input to ChimPipe are **mandatory arguments**. Please see bellow a descripion of them: 

.. code-block:: bash

	-i|--input reads_1.fastq – First mate sequencing reads. ChimPipe deals with paired-end data. 
				   Please make sure the second mate file is in the same directory as 
				   the first one, and the files are named according to the same convention. 
				   E.g: the second mate of "reads_1.fastq" should be "reads_2.fastq". 
						   
	-g|--genome-index genome.gem – Index for the reference genome in GEM format.

	-a|--annotation annotation.gtf – Reference genome annotation file in GTF format. The transcriptome 
						index has to be in the same directory as the annotation. 
								 
	-q|--quality 33 – Quality offset of the FASTQ files [33 | 64 | ignore].

	-l|--seq-library NONE – Type of sequencing library [MATE1_SENSE | MATE2_SENSE | UNSTRANDED]. 
				UNSTRANDED for not strand-specific protocol (unstranded data) and the others for 
				the different types of strand-specific protocols (stranded data).
		          
	-e|--sample-id sample1 – Sample identifier (the output files will be named according to this id).

ChimPipe pipeline has some **optional arguments**, please do ``ChimPipe.sh --help`` to see the list.

.. tip:: If your machine has more than one CPU it is recommended to run ChimPipe with multiple threads. It will speed up the mapping steps a lot. Use the option ``-t|--threads <threads>``, where **threads** is the number of CPUs available. 

.. tip:: It is strongly advisable to use the option ``--similarity-gene-pairs <TEXT FILE>`` to discard junctions connecting genes which encode transcripts with a high sequence homology, which are likely sequencing or mapping artefacts. Please check **Similarity between gene pairs** above to learn how to produce the text file needed. 

.. note:: The pipeline is restartable. That means if ChimPipe fails at some point and you run it again, it will skip the already completed steps. You just need to make sure you remove the files generated in the step the pipeline failed. 

Output
======

By default, ChimPipe produces 3 files as output:

* First mapping BAM file
* Second mapping MAP file
* Chimeric junctions text file

.. tip:: If you want to keep intermediate files in your output run ChimPipe with the flag ``--no-cleanup``. 

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

chr1_121115975_+:chr1_206566046_+	1	1	121115953	206566073	1	1	85450071	GC	AG	SRGAP2D,	SRGAP2,	SRGAP2D,	SRGAP2	.	.	1-1:2,	99.44	1067

