.. _FAQ:

====
FAQ 
====

.. _faq-offset:

Does ChimPipe considers the reads quality for the mapping step? 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The quality score (Q) measures the probability that a base is called incorrectly by the sequencing machine. Within your FASTQ files, they are represented in the fourth line of each read as an ASCII character string (each character corresponds to the Q score of a certain base in the sequencing read). The correspondence between each ASCII character and the Q score is based on some offset. This offset varies depending on the sequencing platform (Illumina machines from CASAVA v1.8 uses 33, while older ones use 64). 

Yes, ChimPipe **deals with both 33 and 64** quality encodings. It will automatically detect it from your reads, so you do not need to specify it as a parameter. 


.. _faq-library:

Which sequencing library protocols are supported?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Different protocols can be used to generate a RNA-seq library. Currently, ChimPipe can handle data generated with the following protocols:

* Non strand-specific protocols. (unstranded data). The information about from which strand the transcript is transcribed is not available. 

* Strand-specific protocols (stranded data):
 
	* **MATE1_SENSE**. Mates 1 are sequenced from the transcript sequence (they will map on the same strand as the transcript), and mates 2 are sequenced from the reverse complement of the transcript sequence (they will map on the strand that is the opposite of the transcript strand). 
	* **MATE2_SENSE**. Mates 1 are sequenced from the reverse complement of the transcript sequence (they will map on the strand that is the opposite of the transcript strand), and mates 2 are sequenced from the transcript sequence (they will map on the same strand as the transcript). 
	
ChimPipe is able to infer the protocol used to produce your data. To do it, it takes a subset of 1M of mapped reads and compares the mapping strand with the strand of the annotated gene where they map. OptionaLly, you can supply this information and skip this step with the option ``-l|--seq-library <library>``.

.. _faq-reads:

Read identifier format
~~~~~~~~~~~~~~~~~~~~~~~

The FASTQ file uses four lines per sequencing read. You need to check the format of the first line of each read, which begins with the '@' character and is followed by a read identifier. This identifier should meet one of the two Illumina standards to specify which member of the pair the read is:

* **CASAVA lower than v1.8**. The identifier has to be a single string ended in /1 or /2 for mate 1 and mate 2, respectively. E. g.:

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
	

* **CASAVA v1.8 or higher**. The identifier consists on two strings separated by a space with the second one specifying the mate in the first digit. E. g:   

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

.. _faq-dependencies:

How can I export the path to the dependencies?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To export the path of bedtools, samtools and blast (if needed) binaries you just need to type:

.. code-block:: bash

	
	$ export PATH=<BEDTOOLS_BINARIES_PATH>:<SAMTOOLS_BINARIES_PATH><BLAST_BINARIES_PATH>:$PATH
	$ # E.g. export bedtools and samtools on my system
	$ export PATH=/users/rg/brodriguez/bin/bedtools2-2.20.1/bin:/users/rg/brodriguez/bin/samtools-0.1.19:$PATH
		

.. _faq-similarity:
How does the script to compute gene pair similarity work?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This script will produce a matrix containing gene pair similarity information through 4 steps:

1. Extract the cDNA sequence of each transcript in the annotation.

2. Make a BLAST database out of the transcript sequences. 

3. Run BLAST on all trancript against all transcripts to detect local similarity between transcripts.

4. Produce a 8 fields matrix where each row corresponds to a gene pair and it contains information about the alignment between the pair of transcripts of this two genes with the maximum alignment similarity and length. Here is a brief description of the 8 fields:

	1. Gene id A
	2. Gene id B
	3. Transcripts alignment similarity
	4. Transcript alignment length
	5. Transcript name A
	6. Transcript name B
	7. Trancript A exonic length
	8. Transcript B exonic length

Note that is expect BLAST binaries to be in your PATH.  

**Example** 

ENSG00000000003.10 ENSG00000003402.15 91.43 70 ENST00000373020.4 ENST00000309955.3 2206 14672


