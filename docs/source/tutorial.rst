.. _tutorial:

=========
Tutorial
=========


The follow tutorial explains how to execute ChimPipe to identify already known fusion genes on the MCF-7 Breast cancer cell line. 

A total of 6 different fusion genes have been experimentally validated by Edgren et al. (Genome Biology, 2011) and Kangaspeska et al. (PLOSone, 2012). 5 of them are detected by ChimPipe:

=============== =================================== ============
 Fusion Gene		Chimeric splice-junction		  Detected
===============	=================================== ============
BCAS4-BCAS3     chr20_49411710_+\:chr17_59445688_+       Yes
ARFGEF2-SULF2   chr20_47538547_+\:chr20_46365686_-       Yes
RPS6KB1-VMP1    chr17_57992064_+\:chr17_57917129_+       Yes
GCN1L1-MSI1     chr12_120628101_-\:chr12_120785317_-     No
AC099850.1-VMP1 chr17_57184952_+\:chr17_57915656_+       Yes
SMARCA4-CARM1   chr19_11097269_+\:chr19_11015627_+       Yes
=============== =================================== ============


1. Download ChimPipe
=====================
Check our :ref:`installation` section for details about how to download and install ChimPipe.

2. Download all-in-one package
===============================

This package contains all the needed files for executing ChimPipe on MCF-7 RNA-seq data in the "input" folder:

* Paired-end RNA-seq data: MCF-7_1.fastq.gz and MCF-7_2.fastq.gz
* Gencode 19 annotation: gencode.v19.annotation.long.gtf
* Genome index: Homo_sapiens.GRCh37.chromosomes.chr.M.gem
* Transcriptome index: gencode.v19.annotation.long.gtf.junctions.gem  
* Transcriptome keys: gencode.v19.annotation.long.gtf.junctions.keys
* Similarity matrix: gencode.v19.annotation.long.similarity.txt

Click `here`_ to download.

.. _here: http://public-docs.crg.es/rguigo/Papers/ChimPipe/ChimPipe_tutorial.tar.gz


3. Execute ChimPipe on MCF-7
=============================

.. code-block:: bash
	
	ChimPipe.sh --fastq_1 MCF-7_1.fastq.gz --fastq_2 MCF-7_2.fastq.gz 
        -g Homo_sapiens.GRCh37.chromosomes.chr.M.gem -a gencode.v19.annotation.long.gtf 
        -t gencode.v19.annotation.long.gtf.junctions.gem -k gencode.v19.annotation.long.gtf.junctions.keys 
        --sample-id MCF-7 --threads 4 --similarity-gene-pairs gencode.v19.annotation.long.similarity.txt
	
4. Examine your output
=======================

ChimPipe identifies 13 chimeric transcript, including 5/6 validated cases and 3 readthrough transcripts. 

We provide ChimPipe output for MCF-7 on the "output" folder. If everything is ok, your results should be identical as the ones provided. 

Check our :ref:`manual` section for an explanation about ChimPipe output files. 



