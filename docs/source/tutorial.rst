.. _tutorial:

=========
Tutorial
=========


The follow tutorial explains how to execute ChimPipe to identify already known fusion genes on the MCF-7 Breast cancer cell line. 

A total of 6 different fusion genes have been experimentally validated by Edgren et al. (Genome Biology, 2011) and Kangaspeska et al. (PLOSone, 2012). 5 of them are detected by ChimPipe:

=============== =================================== ============
 Fusion Gene		Chimeric splice-junction		  Detected
===============	=================================== ============
BCAS4-BCAS3     chr20\_49411710\_+:chr17_59445688\_+       Yes
ARFGEF2-SULF2   chr20\_47538547\_+:chr20_46365686\_-       Yes
RPS6KB1-VMP1    chr17\_57992064\_+:chr17_57917129\_+       Yes
GCN1L1-MSI1     chr12\_120628101\_-:chr12_120785317\_-     No
AC099850.1-VMP1 chr17\_57184952\_+:chr17_57915656\_+       Yes
SMARCA4-CARM1   chr19\_11097269\_+:chr19_11015627\_+       Yes
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



