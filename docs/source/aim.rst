.. _aim.rst:

===============================
ChimPipe for chimeras detection 
===============================

Biological importance of chimeras
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Chimeras are transcripts whose sequence is encoded in two or more different genes. The study of these kind of transcripts is relevant in two different context:

* **Cancer genomics**. It is very well know that the generation of **fusion genes** through chromosomal rearrangements is a major driver in certain types of cancer. These are hydrid genes formed from two previously separate genes that encode altered proteins with abnormal activity. Thus, the identification and study of chimeric transcripts in cancer can serve to detect such alterations in cancer and in last term to the detection of new diagnostic and therapeutic targets. 
  
* **Genome biology**. Understanding how the information is encoded in the genome is one of the challenges of current biology. In the last years, an increasingly number of chimeric transcripts have been reported in non cancer samples from different species. Although their function and the mechanism through they arise are not well know, their existence could reveal new dimensions in how information is encoded in the genome. As, the combination of sequences from different genes could be a way of increasing the information content of the genomes. Several **transcriptional events** have been proposed to explain how this transcripts could arise:

	* **Trans-splicing**. Splicing between two pre-mRNAs from different genes producing a chimera. 
	
	* **Readthrough or conjoined transcription**. Combined transcription of two consecutive genes on the same chromosome and strand. Typically, such chimera begin at the promoter of the upstream gene and end at the termination point of the downstream gene.
	
	* **SHS-transcriptional slippage**. Pairing of short homologous sequences (SHS) among two different transcripts that leads to a recombination producing a chimeric transcript. 


ChimPipe principle
~~~~~~~~~~~~~~~~~~~
Briefly, ChimPipe uses the `GEMtools RNA-seq pipeline`_ to continuously map the read pairs to the genome and the annotated transcriptome, and segmentally maps them into the genome to find *de novo* splice junctions in the same chromosome and strand. Then, in a second mapping step, it further segmentally maps the reads that could not be mapped this way to find *de novo* splice junctions allowing for different chromosome and strand through the `GEM RNA mapper`_. ChimPipe further aggregates the two block of segmentally mapped reads whose two parts map in two different genes into chimeric junctions, and applies several filters such as number of reads spanning and read pairs encompassing the chimeric junction. 

Currently, we are working in ChimPipe's paper where our approach is deeply described. 

.. _GEMtools RNA-seq pipeline: http://gemtools.github.io/
.. _GEM RNA mapper: http://algorithms.cnag.cat/wiki/The_GEM_library



ChimPipe features
~~~~~~~~~~~~~~~~~~

* Can handle both unstranded and stranded data
* Provides the precise chimeric junction coordinates
* Outputs many pieces of information for each chimeric junction
* Has a good balance sensitivity/specificity according to our benchmark



