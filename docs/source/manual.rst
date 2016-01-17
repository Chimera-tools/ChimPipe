.. _manual:

======
Manual
======

This section describes ChimPipe' input files, how to generate them, how to run ChimPipe, and how to interpret its outputs. 

Input files
===========
ChimPipe takes as input 4 mandatory files and 1 optional file.  

Mandatory:

* :ref:`Paired-end (PE) RNA-seq reads (FASTQ or BAM) <input-reads>`
* :ref:`Genome index <input-genome-index>` 
* :ref:`Reference annotation <input-annot>`
* :ref:`Transcriptome index <input-annot-index>`

Optional:


* :ref:`Gene pair similarity file <input-similarity>`

.. _input-reads:

Paired-end (PE) RNA-seq reads (FASTQ or BAM)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ChimPipe has been designed to deal with `Illumina paired-end`_ RNA sequencing data. It can take as input 2 `FASTQ`_ files (one for each member in the pair) or a single `BAM`_ file with already mapped reads.

.. _Illumina paired-end: http://technology.illumina.com/technology/next-generation-sequencing/paired-end-sequencing_assay.ilmn
.. _FASTQ: http://maq.sourceforge.net/fastq.shtml
.. _BAM: http://samtools.github.io/hts-specs/SAMv1.pdf

.. warning:: **FASTQ as input.** The read identifier should follow Illumina convention to specify which member of the pair the read is. See our :ref:`FAQ <faq-reads>` section for more information. 

.. warning:: **BAM as input.** Your BAM file should contain the **NH** (Number of hits in the reference) and **NM** (Edit distance to the reference) optional fields. Most popular mappers produce a BAM with these fields. Please check `SAM`_ specifications for more information. 

.. _SAM: http://samtools.github.io/hts-specs/SAMv1.pdf

.. _input-genome-index:

Genome index
~~~~~~~~~~~~
An indexed reference genome in GEM format has to be provided to do the mapping steps. 

We provide some **pre-generated genome indices** for human, mouse and drosophila in the :ref:`Downloads` section, however if you have a different genome, you just need to run the *GEMtools indexer* (at ``ChimPipe/bin/gemtools-1.7.1-i3/gemtools``) with your genome in FASTA format to produce it:

.. _FASTA:
 
.. code-block:: bash

	$ gemtools index -i genome.fa 

It will produce 3 files in the directory where the genome fasta file is located:

* **genome.gem** – indexed genome file in GEM format (needed for running ChimPipe).   
* genome.hash – hash table with the genome (not needed). 
* genome.log – indexer log file (not needed).    

.. tip:: If your machine has more than one CPU, it is recommended to run the indexer with multiple threads. Use the option ``-t <threads>``, where **threads** is the number of available CPUs. 

.. _input-annot:

Reference annotation
~~~~~~~~~~~~~~~~~~~~~
ChimPipe also takes as input a reference annotation in `GTF`_ format (it has to include annotated exons). It can contain other features different from exons, i. e. introns or UTR, but they will not be considered by ChimPipe in the chimera detection process. This annotation needs to contain at least two (tag,value) pairs in the attribute field with the "gene_id" and the "transcript_id" tags. Two optional (tag,value) pairs will be taken into account by ChimPipe if they are provided: "gene_name" and "gene_type". e.g:

.. _GTF: http://www.ensembl.org/info/website/upload/gff.html

.. code-block:: bash
	
	# This is an example of an annotated exon with an appropiate format. 	
	# The attributes are the gene_id, transcript_id (mandatory), the gene type and gene name (optional), 
	# plus some additional (tag,value) pairs that will not be considered by ChimPipe.   
	
	chr1	HAVANA	exon	69091	70008	.	+	.	gene_id "ENSG00000186092.4"; transcript_id "ENST00000335137.3"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "OR4F5";
	transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "OR4F5-001"; exon_number 1; exon_id "ENSE00002319515.1"; level 2; tag "basic"; tag "appris_principal";
	tag	"CCDS"; ccdsid "CCDS30547.1"; havana_gene "OTTHUMG00000001094.1"; havana_transcript "OTTHUMT00000003223.1";

.. note:: ChimPipe has been benchmarked with `Gencode v10`_ and `UCSC Known Genes`_  human gene annotations. It displayed a better sensitivity with Gencode v10 but a similar false positive rate with both annotations. It is advisable to use Gencode annotation, since it is a richer annotation which increases the sensitivity of the chimera detection process. 

.. _Gencode v10: http://www.gencodegenes.org/releases/10.html
.. _UCSC Known Genes: https://genome.ucsc.edu/cgi-bin/hgTables?command=start

.. _input-annot-index:

Transcriptome index
~~~~~~~~~~~~~~~~~~~
A transcriptome index in GEM format has to be provided as input to ChimPipe to find reads spanning annotated splice junctions. 

We provide some **pre-generated transcriptome indices** for human, mouse and drosophila annotations in the :ref:`Downloads` section, however if your genome annotation or your genome is different, you will need to to run the *GEMtools transcriptome indexer* ((at ``ChimPipe/bin/gemtools-1.7.1-i3/gemtools``)) on your previously generated GEM indexed genome and your annotation in GTF format, as indicated below. 

.. code-block:: bash

	$ gemtools t-index -i genome.gem -a annotation.gtf	

It will produce 5 files in your current working directory:

* annotation.gtf.junctions – annotated splice junctions coordinates (not needed)
* annotation.gtf.junctions.fa – annotated splice junctions sequence (not needed)
* **annotation.gtf.junctions.gem** – transcriptome index in GEM format (needed)
* **annotation.gtf.junctions.keys** – keys to convert from transcriptome to genome (needed)
* annotation.gtf.junctions.log – indexer log file (not needed)

.. tip:: If your machine has more than one CPU it is recommended to run the indexer with multiple threads. Use the option ``-t <threads>``, where **threads** is the number of available CPUs. 

.. _input-similarity:

Gene pair similarity file (Optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
One of ChimPipe's filtering steps to discard actefactual chimeras is to filter out those chimeric junctions connecting genes that encode transcripts with a high sequence homology. To do this ChimPipe produces a matrix with information about the sequence similarity between annotated gene pairs. This step is done internally and takes around 1:30h. So, in case you want to run many samples with the same annotation, it is recommended to provide the matrix with the option ``--similarity-gene-pairs <MATRIX TEXT FILE>``, to avoid recomputing it for every sample. 

You can download our **pre-generated matrices** por human, mouse and drosophila annotations from :ref:`Downloads` section or you can produce your own matrix with the script ``ChimPipe/src/bash/similarity_bt_gnpairs.sh`` as follows:

        $ bash similarity_bt_gnpairs.sh annot.gtf genome.gem

Please check out our :ref:`FAQ <faq-similarity>` section for more information about how the script works. 

.. warning:: Make sure you run ChimPipe with a similarity matrix generated from the annotation and genome you are using.  


Executing ChimPipe
==================

1. Setting up the environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
As explained in the :ref:`installation` section, you need to have BEDtools, SAMtools and Blast installed on your system to execute ChimPipe. In case you do not have them, you can download and install them from their web pages. 

Once they are installed, you can export the path to their binaries with a simple bash command:

.. code-block:: bash

	
	$ export PATH=<BEDTOOLS_BINARIES_PATH>:<SAMTOOLS_BINARIES_PATH><BLAST_BINARIES_PATH>:$PATH
	$ # E.g. export bedtools and samtools on my system
	$ export PATH=~/bin/bedtools2-2.20.1/bin:~/bin/samtools-0.1.19:$PATH


2. Running ChimPipe
~~~~~~~~~~~~~~~~~~~
There are two different running modes depending if what you will provide FASTQ or BAM files:

A) FASTQ
---------

.. code-block:: bash
	
	ChimPipe.sh --fastq_1 <mate1_fastq> --fastq_2 <mate2_fastq> -g <genome_index> -a <annotation> -t <transcriptome_index> -k <transcriptome_keys> [OPTIONS]

All these files and parameters given as input to ChimPipe are **mandatory arguments**. Please see bellow their descripion: 

.. code-block:: bash

        --fastq_1                       <FASTQ>         First mate sequencing reads in FASTQ format. It can be gzip compressed [.gz].
        --fastq_2                       <FASTQ>         Second mate sequencing reads in FASTQ format. It can be gzip compressed [.gz].
        -g|--genome-index               <GEM>           Reference genome index in GEM format.
        -a|--annotation                 <GTF>           Reference gene annotation file in GTF format.                                
        -t|--transcriptome-index        <GEM>           Annotated transcriptome index in GEM format.
        -k|--transcriptome-keys         <KEYS>          Transcriptome to genome indices coordinate conversion keys (generated when producing transcriptome index).  
        --sample-id                     <STRING>        Sample identifier (output files are named according to this id).  

B) BAM
--------

.. code-block:: bash

        ChimPipe.sh --bam <bam> -g <genome_index> -a <annotation> [OPTIONS]

.. code-block:: bash

        --bam                           <BAM>           Already mapped reads in BAM format, then first mapping step skipped.
        -g|--genome-index               <GEM>           Reference genome index in GEM format.
        -a|--annotation                 <GTF>           Reference genome annotation file in GTF format.
        --sample-id                     <STRING>        Sample identifier (the output files will be named according to this id).  

.. warning:: ChimPipe results are highly dependent on the mapping step. BAM files should have been generated with a mapper able to align read-pairs across different chromosomes and strands. If you are unsure about that, please convert your BAM file to FASTQ and run ChimPipe in FASTQ as input mode. 
								 
**Optional arguments.** Do ``ChimPipe.sh -h or --help`` to see a short help with the main options. You can also do ``ChimPipe.sh --full-help`` to see the all the possible 
options. 

.. tip:: If your machine has more than one CPU it is recommended to **run ChimPipe with multiple threads** (at least 4). It will speed up the mapping steps a lot. Use the option ``-t|--threads <threads>``, where threads is the number of CPUs available. 

.. note:: The **pipeline is restartable**: this means that if ChimPipe fails and you run it again, it will skip the steps that have been already completed. You just need to make sure that you removed the files generated in the step where the pipeline failed. 

Output
======

By default, ChimPipe produces 4 main output files:

* :ref:`First mapping BAM <output-bam>` (MappingPhase/FirstMapping/[sample_id]_firstMap.bam).
* :ref:`Second mapping MAP <output-map>` (MappingPhase/SecondMapping/[sample_id]_secondMap.map).
* :ref:`Final chimeric junctions <output-chimeras>` (chimericJunctions_[sample_id].txt).
* :ref:`Discarded chimeric junctions <output-chimeras>` (chimericJunctions_filtered_[sample_id].txt).

.. tip:: If you want to keep intermediate output files, run ChimPipe with the ``--no-cleanup`` option. 

First mapping BAM file
~~~~~~~~~~~~~~~~~~~~~~
`BAM`_ file containing the reads mapped in the genome, transcriptome and *de novo* transcriptome with the `GEMtools RNA-seq pipeline`_. 

This is the standard format for next-generation sequencing, meaning that most analysis tools work with this format. The bam file produced can therefore be used to do other downstream analyses such as gene and transcript expression quantification.

.. _BAM: http://samtools.github.io/hts-specs/SAMv1.pdf
.. _GEMtools RNA-seq pipeline: http://gemtools.github.io/

Second mapping MAP file
~~~~~~~~~~~~~~~~~~~~~~~
MAP file containing reads segmentally mapped in the genome allowing for interchromosomal, different strand and unexpected genomic order mappings. 

Final and filtered chimeric junction files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Two tabular text files with the detected and discarded chimeric splice junctions from your RNA-seq dataset. They consist on rows of 35 fields, where each row corresponds to a chimeric junction and each field contains a piece of information about the chimera. Here is a brief description of the 35 fields (most relevant fields highlighted in bold):

1. **juncCoord** - Position of the chimeric splice junction in the genome described as follows: chrA"_"coordA"_"strandA":"chrB"_"coordB"_"strandB. E. g., "chr4_90653092_+:chr17_22023757_-" is a chimeric junction between the position 90653092 of chromosome 4 in plus strand, and the position 22023757 of chromosome chr17 in minus strand. Junction coordinates defined using 1-based system.

2. **type** - Type of chimeric splice junction. Junctions classified in 5 different categories:
	
	- readthrough: donor and acceptor sites within the same chromosome, strand and within less than 100.000 base pairs.
	- intrachromosomal: donor and acceptor sites within the same chromosome, strand and in separated by more than 100.000 base pairs.
	- inverted: donor site downstream than acceptor site (opposite as expected) and both in the same chromosome.
	- interstrand: donor and acceptor sites within the same chromosome but in different strand.
	- interchromosomal: donor and acceptor sites in different chromosomes.
	
3. **filtered** - Flag to specify if the chimeric junction has been filtered out (1) or not (0).

4. **reason** - List of filters the chimeric junction failed to pass. There is a tag per filter:
	
	- totalSupport.
	- spanningReads.
	- consistentPE.
	- totalSupport.
	- spanningReads.
	- consistentPE.
	- percStag.
	- percMulti.
	- percInconsistentPE.
	- similarity.
	- biotype.
	
5. **nbTotal(spanning+consistent)** - Total supporting evidences (spanning reads + consistent paired-ends).

6. **nbSpanningReads** - Number of split-reads spanning the chimeric splice junction.  

7. nbStaggered - Number of spanning split-reads aligning at different positions.

8. percStaggered - Percentage of spanning split-reads that aligns at different positions.  

9. nbMulti - Number of multimapped split-reads spanning the chimeric junction. 

10. percMulti - Percentage of spanning split-reads mapping in multiple locations. 

11. **nbConsistentPE** - Number of discordant paired-end consistent with the chimeric junction.

12. nbInconsistentPE - Number of discordant paired-end inconsistent with the chimeric junction.

13. percInconsistentPE - Percentage of discordant paired-end inconsistent with the chimeric junction. 

14. overlapA - Percentage of overlap between the 5' split-read cluster and the annotated exon.

15. overlapB - Percentage of overlap between the 3' split-read cluster and the annotated exon

16. **distExonBoundaryA** - Distance between the chimeric junction donor site and the exon boundary (annotated donor). 

17. **distExonBoundaryB** - Distance between the chimeric junction acceptor site and the exon boundary (annotated acceptor).

18. blastAlignLen - Maximum length of the BLAST alignment between all the transcripts of the gene pairs connected by the chimeric junction. ”na” if no blast hit found. 

19. blastAlignSim - Maximum percent of similarity in the BLAST alignment between the transcript with the longest BLAST alignment. ”na” if no blast hit found.

20. **donorSS**	- Splice donor site sequence.

21. **acceptorSS** - Slice acceptor site sequence.

22.	beg	- Split-reads cluster 5' coordinates. 

23. end	- Split-reads cluster 3' coordinates.

24. sameChrStr - Flag to specify if the connected gene pairs are in the same cromosome and strand (1) or not (0).	

25. okGxOrder -	Flag to specify if the connected gene pairs are in genomic order (1) or not (0). "na" in case the samechrstr field was 0 (being in genomic order means that the donor gene is located in 5' while the acceptor in 3'. This means that if both genes are in the plus strand, the genomic coordinates of the first gene are lower than the ones for the second one. For genes in the minus it is the opposite).

26. dist -	Distance between the chimeric junction splice sites. "na" in case the “sameChrStr” field was 0.

27. **gnIdsA** - 5' gene/s involved in the chimeric transcript. 	

28. **gnIdsB** - 3' gene/s involved in the chimeric transcript. 		

29. **gnNamesA** - Name of the genes in the field *gnIdsA*. "na" if unknown.

30. **gnNamesB** - Name of the genes in the field *gnIdsB*. "na" if unknown.	

31. **gnTypesA** - Biotype of the genes in the field *gnIdsA*. "na" if unknown. 	

32. **gnTypesB** - Biotype of the genes in the field *gnIdsB*. "na" if unknown. 		

33. juncSpanningReadsIds - Identifiers of the split-reads spanning the chimeric splice junction. 

34. consistentPEIds	- Identifiers of the paired-ends consistent with  the chimeric splice junction. 

35. inconsistentPEIds - Identifiers of the paired-ends inconsistent with  the chimeric splice junction. 


