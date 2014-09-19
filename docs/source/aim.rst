.. _aim.rst:

===============================
ChimPipe for chimeras detection 
===============================

Biological importance of chimeras
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Chimeras are transcripts whose sequence is encoded in two or more different genes. The study of these kind of transcripts is relevant in two different context:

* **Cancer genomics**. It is very well know that the generation of **fusion genes** through chromosomic rearrangements is a major driver in certain types of cancer. These are hydrid genes formed from two previously separate genes that encode altered proteins with abnormal activity. Thus, the identification and study of chimeric transcripts in cancer can serve to detect such alterations in cancer and in last term to the detection of new diagnostic and therapeutic targets. 
  
* **Genome biology**. Understanding how the information is encoded in the genome is one of the challenges of current biology. In the last years, an increasingly number of chimeric transcripts have been reported in non cancer samples from different species. Although their function and the mechanism through they arise are not well know, their existence could reveal new dimensions in how information is encoded in the genome. As, the combination of sequences from different genes could be a way of increasing the information content of the genomes. Several **transcriptional events** have been proposed to explain how this transcripts could arise:

	* **Trans-splicing**. Splicing between two pre-mRNAs from different genes producing a chimera. 
	
	* **Readthrough or conjoined transcription**. Combined transcription of two consecutive genes on the same chromosome and strand. Typically, such chimera begin at the promoter of the upstream gene and end at the termination point of the downstream gene.
	
	* **SHS-transcriptional slippage**. Pairing of short homologous sequences (SHS) among two different transcripts that leads to a recombination producing a chimeric transcript. 


ChimPipe features
~~~~~~~~~~~~~~~~~~

* Aims to detect chimeric transcripts arising from transcriptional events and genomic rearrangements.

* Deals with stranded and unstranded **Paired End RNA-seq** data. 	

* Provides precise chimeric junction coordinates between connected genes.

* Ouput many qualitative and quantitative information associated to the chimeric junctions.

* Considers two sources of evidence to find them: reads spanning and read pairs encompassing chimeric junctions.

* Performs a continuous and segmental aligment of the reads to detect them. 

* Read mapping done with `GEM`_ and the `GEM rna-mapper`_, two exhaustive read mappers, through the `GEMtools`_ library. 

* Implements a set of filters to discard spurious chimeras arising from sequencing and alignment artefacts: 

	* Nb. staggered spanning reads filter
	* Nb. encompassing read pairsv filter
	* Sequence similarity filter between connected gene pairs
	* PCR-duplicates filter
	* Minimum anchor length filter

* Evaluated on several positive and negative datasets on which other programs were already benchmarked, and out of five programs, ChimPipe ranked first for sensitivity, third for specificity and first overall. Unpublished data, we are working on the manuscript.  

.. _GEM: http://algorithms.cnag.cat/wiki/The_GEM_library
.. _GEM rna-mapper: http://algorithms.cnag.cat/wiki/The_GEM_library
.. _GEMtools: http://gemtools.github.io/

