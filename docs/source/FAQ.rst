.. _FAQ:

====
FAQ 
====

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

.. _faq-similarity:

How the script to compute gene pair similarity works?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

**Example** 

ENSG00000000003.10 ENSG00000003402.15 91.43 70 ENST00000373020.4 ENST00000309955.3 2206 14672

.. _faq-dependencies:

How can I export the path to the dependencies?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To export the path of bedtools, samtools and blast (if needed) binaries you just need to type:

.. code-block:: bash

	
	$ export PATH=<BEDTOOLS_BINARIES_PATH>:<SAMTOOLS_BINARIES_PATH><BLAST_BINARIES_PATH>:$PATH
	$ # E.g. export bedtools and samtools on my system
	$ export PATH=/users/rg/brodriguez/bin/bedtools2-2.20.1/bin:/users/rg/brodriguez/bin/samtools-0.1.19:$PATH
	
.. _faq-offset:

How can I know the quality offset of my RNA-seq reads?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We provide a bash script to detect the offset of your RNA-seq data. You can find it at ``ChimPipe/tools/detect.fq.qual.sh``.

.. code-block:: bash

	$ bash detect.fq.qual.sh reads_1.fastq
	$ Offset 33

.. _faq-library:
	
How can I know the RNA-seq library type?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We supply a bash script to infer the protocol used to produce the sequencing library from a subset of mapped read in BAM format. It is at ``ChimPipe/tools/infer_exp.sh``. You can proceed as follow:

.. code-block:: bash

	$ # 1) First you need to extract and map a subset of read pairs (E.g 1000 pairs)
	$ head -4000 reads_1.fastq > reads.subset_1.fastq
	$ head -4000 reads_2.fastq > reads.subset_2.fastq
	$ # 2) Map the reads with the gemtools rna-pipeline (at `ChimPipe/bin/gemtools-1.7.1-i3/gemtools`)	
	$ gemtools rna-pipeline -f reads.subset_1.fastq reads.subset_2.fastq -i genome.gem -a annotation.gtf -q 33 
	$ # Note: you also need your transcriptome index in the same directory as the annotation  
	$ # 3) Run our script with the BAM file generated and the annotation used for the mapping step.
	$ bash infer_exp.sh annotation.gtf reads_subset.bam
	$ bash infer_exp.sh annotation.gtf reads_subset.bam

This script will produce an output like this:

.. code-block:: bash

	Creating the reference gene model in bed format

	Inferring experiment
	Reading reference gene model gen10.long.bed ... Done
	Loading SAM/BAM file ...  Total 200000 usable reads were sampled

	This is PairEnd Data
	Fraction of reads explained by "1++,1--,2+-,2-+": 0.0135 
	Fraction of reads explained by "1+-,1-+,2++,2--": 0.9866  
	Fraction of reads explained by other combinations: 0.0000
	Removing temporary files
	DONE

The three rows starting with "Fraction of reads explained by.." are the ones giving information regarding the library. They contain several strings of three characters, i.e. 1+-, where:

* Character 1. 1 and 2 are mate1 and mate2 respectively. 
* Character 2. + and - is the strand where the read maps. 
* Character 3. + and - is the strand where the gene in which the read overlaps is annotated. 

You can apply the following rules to infer the used library from this information:

* **NONE**. Not strand-specific protocol (unstranded data). Fraction of reads explained by "1++,1--,2+-,2-+" and "1+-,1-+,2++,2--" close to 0.5000 in both cases. 

Strand-specific protocols (stranded data):
 
* **MATE1_SENSE**. Fraction of reads explained by "1++,1--,2+-,2-+" close to 1.0000. 
* **MATE2_SENSE**. Fraction of reads explained by "1+-,1-+,2++,2--" close to 1.0000.

.. note:: This script makes use of a tool from the `RSeQC package`_. So, you should have RSeQC installed in your system to run the script. We plan to recode it soon to use BEDtools instead of RSeQC and to avoid having an additional dependency.  

.. _RSeQC package: http://rseqc.sourceforge.net/.

.. tip:: In case you have any problem to map the subset of reads you can check the `GEMtools RNA Pipeline Quickstart`_ for more details. 

.. _GEMtools RNA Pipeline Quickstart: http://gemtools.github.io/docs/rna_pipeline.html


