.. _FAQ:

====
FAQ 
====

.. _faq-offset:

Quality offset encoding 
~~~~~~~~~~~~~~~~~~~~~~~~
ChimPipe accepts FASTQ with both 33 and 64 (old Illumina machines) quality encodings. It will automatically infer the offset from your input reads, so you do not need to specify it as a parameter. 

.. _faq-library:

RNA-seq protocols accepted
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Different protocols can be used to generate a RNA-seq library. Currently, ChimPipe can handle data generated with the following protocols:

* **Non strand-specific protocols**. Information about from which strand the transcript is transcribed is not available. 

* **Strand-specific protocols**:
 
	* **MATE1_SENSE**. Mate 1 directly sequenced from the transcript sequence and mate 2 from the reverse complement of the transcript sequence. 
	* **MATE2_SENSE**. Mate 1 sequenced from the reverse complement of the transcript sequence and mate 2 directly from the transcript sequence. 
	
ChimPipe is able to infer the protocol used to produce your data. It takes a subset of 1M of mapped reads and compares the mapping strand with the strand of the annotated gene where they map. Optionally, you can supply this information and skip this step with the option ``-l|--seq-library <library>``.

.. _faq-reads:

Read identifier format
~~~~~~~~~~~~~~~~~~~~~~~

ChimPipe requires the read identifiers to meet one of the two Illumina standards to specify which member of the pair the read is:

* **CASAVA lower than v1.8**. The identifier has to end in /1 and /2 for mate 1 and mate 2, respectively. E. g.:

.. code-block:: bash
	
	$ # Mate 1
	@SRR018259.1_BI:080831_SL-XAN_0004_30BV1AAXX:5:1:708:1868#0/1
	GTAACATATTCACAGACATGTCAGTGTGCACTCAGGAACACTGGTTTCATT
	+
	IIIIIIIIIIIIIII>=CII=8=H032/-++D+'.@)2+4/+)1'4.#"*.
	
	$ # Mate 2
	@SRR018259.1_BI:080831_SL-XAN_0004_30BV1AAXX:5:1:708:1868#0/2
	CAGATGGCATTGCAGAACTAAAAAAGCAAGCTGAATCAGTGTTAGATCTCC
	+
	IIIGIIIIIIIIIFI:>DID<AFH=>78@I41I055549.?+.42-'%**'
	

* **CASAVA v1.8 or higher**. The identifier consists on two strings separated by a space with the second one specifying the mate in the first digit. E.g:   

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

In case your FASTQ files do not meet this format you should modify the identifier. Awk is a perfect tool for such kind of problems. E.g: we downloaded a pair of FASTQ files from the NCBI Sequence Read archive:

.. code-block:: bash
	
	# Mate 1 reads without a proper identifier.  
	
	@SRR018259.1 BI:080831_SL-XAN_0004_30BV1AAXX:5:1:708:1868 length=51
	GTAACATATTCACAGACATGTCAGTGTGCACTCAGGAACACTGGTTTCATT
	+SRR018259.1 BI:080831_SL-XAN_0004_30BV1AAXX:5:1:708:1868 length=51
	IIIIIIIIIIIIIII>=CII=8=H032/-++D+'.@)2+4/+)1'4.#"*.

	# No worries, we can use awk to fix the FASTQ file
	
	$ awk '{if (NR%2==1){print "@"$2"#0/1"} else {print $0}}' dataset_1.fastq 		
	
	@BI:080831_SL-XAN_0004_30BV1AAXX:5:1:708:1868#0/1
	GTAACATATTCACAGACATGTCAGTGTGCACTCAGGAACACTGGTTTCATT
	@BI:080831_SL-XAN_0004_30BV1AAXX:5:1:708:1868#0/1
	IIIIIIIIIIIIIII>=CII=8=H032/-++D+'.@)2+4/+)1'4.#"*.

	$ # Finally, Apply the same procedure for mate 2 FASTQ

.. _faq-dependencies:



