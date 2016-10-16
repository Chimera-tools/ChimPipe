.. _manual:

======
Manual
======

This section describes ChimPipe input files, how to run ChimPipe and interpret the output files. 

Input 
======
ChimPipe takes as input 3-4 mandatory files and 1 optional file.  

Mandatory:

* :ref:`Paired-end RNA-seq reads<input-reads>`
* :ref:`Genome index <input-genome-index>` 
* :ref:`Reference annotation <input-annot>`
* :ref:`Transcriptome index <input-annot-index>` (only in FASTQ as input mode)

Optional:

* :ref:`Gene pair similarity <input-similarity>`

.. _input-reads:

Paired-end RNA-seq reads (FASTQ or BAM)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ChimPipe makes use of `Illumina`_ paired-end RNA sequencing data to identify chimeric transcripts. It can deal with both unaligned and pre-aligned reads in the standard formats:

A) **FASTQ**:
-----------------

* Two independent FASTQ files are required for mate 1 and 2 reads, respectively. 
* The read identifiers should follow any of the 2 `Illumina conventions`_ to specify which member of the pair the reads are. Check our :ref:`FAQ` section for more information.

B) **BAM**:  
---------------

* In BAM as input mode the first mapping step is skipped. 
* A single BAM file with pre-aligned reads has to be given as input. 
* The BAM file should contain unmapped reads. ChimPipe will attempt to split-map them across different chromosomes and strands to identify chimeric splice junctions in a second mapping step. 
* The **NH** (Number of hits) and **NM** (Edit distance) fields are required. Most popular mappers provide this information. Check `SAM`_ specifications for more information. 


For more information about the supported input RNA-seq data check our :ref:`FAQ` section. 

.. note:: ChimPipe results are highly dependent on the mapping step. BAM files should have been generated with a mapper able to align read-pairs across different chromosomes and strands. If you are unsure about that, please convert your BAM file to FASTQ and run ChimPipe in FASTQ as input mode. 
		
.. _Illumina: http://technology.illumina.com/technology/next-generation-sequencing/paired-end-sequencing_assay.ilmn
.. _FASTQ: http://maq.sourceforge.net/fastq.shtml
.. _Illumina conventions: https://en.wikipedia.org/wiki/FASTQ_format
.. _BAM: http://samtools.github.io/hts-specs/SAMv1.pdf
.. _SAM: http://samtools.github.io/hts-specs/SAMv1.pdf

.. _input-genome-index:

Genome index
~~~~~~~~~~~~
An indexed reference genome in GEM format is required for first and second mapping steps.

We supply **pre-generated indices** for Human, Mouse and Drosophila in the :ref:`Downloads` section. 

If you aim to execute ChimPipe with your own genome you just need to use the *GEMtools indexer* (at ``ChimPipe/bin/gemtools-1.7.1-i3/gemtools``):

.. _FASTA:
 
.. code-block:: bash

	$ gemtools index -i genome.fa 

It will produce 3 files in the directory where your genome is located:

* **genome.gem** – indexed genome in GEM format (mandatory).   
* genome.hash – genome hash table. 
* genome.log – log file.    

.. tip:: If your machine has more than one CPU, it is recommended to run the indexer with multiple threads. Use the option ``-t <threads>``, where **threads** is the number of available CPUs. 

.. _input-annot:

Reference annotation
~~~~~~~~~~~~~~~~~~~~~
ChimPipe requires a reference gene annotation in the standard `GTF`_ format:

* The annotation has to contain all the annotated exons for a given species. 
* Any other feature is accepted, Introns, UTR.., but they will not be considered for chimera detection.
* The gtf file has to be sorted by chr, then start and end.
* The following mandatory tag,value pairs are required in the attribute field: gene_id "X"; transcript_id "Y";
* Two optional tag,value pairs will be considered if provided: "gene_name" and "gene_type". There can be additional tag,value pairs, but they will not be taken into account.
* The order of the pairs does not matter, but the tag ids have to be exactly the ones as described since they are used to parse the gtf and extract the information. 

E.g:

.. code-block:: bash
	
	# Example of an annotated human exon in GTF format. 	
	# The attributes are gene_id (mandatory), gene type and gene name (optional) 
	# plus some additional "tag,value" pairs that will not be considered by ChimPipe.   
	
	chr1	HAVANA	exon	69091	70008	.	+	.	
        gene_id "ENSG00000186092.4"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "OR4F5";

.. tip:: We have extensively tested and applied ChimPipe to analyse human and mouse RNA-seq data with `Gencode`_ annotations. So, we suggest to use Gencode if possible.   

.. _GTF: http://www.ensembl.org/info/website/upload/gff.html
.. _Gencode: http://www.gencodegenes.org/

.. _input-annot-index:

Transcriptome index (only in FASTQ as input mode)
~~~~~~~~~~~~~~~~~~~
A transcriptome index containing all the annotated exon-exon junctions for each transcript is needed for the first mapping step. 

We supply **pre-generated transcriptome indices** for Human, Mouse and Drosophila annotations in the :ref:`Downloads` section.

If you aim to execute ChimPipe with your own annotation you just need to use the *GEMtools transcriptome indexer* (at ``ChimPipe/bin/gemtools-1.7.1-i3/gemtools``) on your previously generated GEM indexed genome and annotation, as indicated below:

.. code-block:: bash

	$ gemtools t-index -i genome.gem -a annotation.gtf	

It will produce 5 files in your current working directory:

* annotation.gtf.junctions – annotated splice junctions coordinates. 
* annotation.gtf.junctions.fa – annotated splice junctions sequences. 
* **annotation.gtf.junctions.gem** – transcriptome index in GEM format.
* **annotation.gtf.junctions.keys** – keys to convert from transcriptome to genome coordinates. 
* annotation.gtf.junctions.log – log file. 

.. tip:: If your machine has more than one CPU it is recommended to run the indexer with multiple threads. Use the option ``-t <threads>``, where **threads** is the number of available CPUs. 

.. _input-similarity:

Gene pair similarity file (Optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ChimPipe filters out artefactual chimeric junctions involving genes with high exonic sequence homology. These genes are prone to produce spurious misaligned split-reads connecting them. 

Before applying this filter, ChimPipe has to compute a similarity matrix between every annotated gene pair. 

This step takes around 45'-60' depending on the annotation. So, in case you plan to run many samples with the same annotation, it is recommended to choose one of these two options:

**A)** Execute ChimPipe with a single sample and then reuse the generated matrix (``${outDir}/GnSimilarity/${annotation_id}.similarity.txt``) to run the other samples 

**B)** Pre-compute the matrix executing ``ChimPipe/src/bash/similarity_bt_gnpairs.sh`` as follows and run all your samples:

        $ bash similarity_bt_gnpairs.sh annot.gtf genome.gem

Either if you use A) or B) you can provide the matrix to ChimPipe with the option ``--similarity-gene-pairs <MATRIX TEXT FILE>``. Otherwise, ChimPipe will generate the same matrix per each sample.

Note, we supply **pre-generated matrices** for Human, Mouse and Drosophila in the :ref:`Downloads` section. 

.. warning:: Make sure you run ChimPipe with a similarity matrix generated from the same reference annotation and genome you are using.  


Execute ChimPipe
================

1. Setting up the environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
As explained in the :ref:`installation` section, to execute ChimPipe you need to have BEDtools, SAMtools and Blast binaries installed in your linux environment. 

.. code-block:: bash

	$ # Example about how to export your binaries under your environmnet
	$ export PATH=<BEDTOOLS_PATH>:<SAMTOOLS_PATH><BLAST_PATH>:$PATH
	$ export PATH=~/bin/bedtools2-2.20.1/bin:~/bin/samtools-0.1.19:~/bin/blastn:$PATH

2. Running ChimPipe
~~~~~~~~~~~~~~~~~~~
ChimPipe has two different running modes:

A) FASTQ
---------

.. code-block:: bash
	
	ChimPipe.sh --fastq_1 <mate1_fastq> --fastq_2 <mate2_fastq> -g <genome_index> 
        -a <annotation> -t <transcriptome_index> -k <transcriptome_keys> [OPTIONS]

All these files and parameters given as input to ChimPipe are **mandatory arguments**. See bellow a detailed descripion: 

.. code-block:: bash

        --fastq_1                       <FASTQ>         First mate sequencing reads in FASTQ format. 
                                                        It can be gzip compressed [.gz].
        --fastq_2                       <FASTQ>         Second mate sequencing reads in FASTQ format. 
                                                        It can be gzip compressed [.gz].
        -g|--genome-index               <GEM>           Reference genome index in GEM format.
        -a|--annotation                 <GTF>           Reference gene annotation file in GTF format.                                
        -t|--transcriptome-index        <GEM>           Annotated transcriptome index in GEM format.
        -k|--transcriptome-keys         <KEYS>          Transcriptome to genome indices coordinate conversion keys 
                                                        (generated when producing transcriptome index).  
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

						 
**Optional arguments.** Do ``ChimPipe.sh -h or --help`` to see a short help with the main options. You can also do ``ChimPipe.sh --full-help`` to see the all the possible 
options. 

.. tip:: If your machine has more than one CPU it is recommended to **run ChimPipe with multiple threads** (at least 4). It will speed up the mapping steps. Use the option ``-t|--threads <threads>``, where threads is the number of CPUs available. 

.. note:: **Checkpoints**: ChimPipe has a checkpoint in every step. So, in case the program fail at one step. You can restart the analysis from this step. You just need to make sure that you removed the files generated in the step where the pipeline failed. 

Output
======

By default, ChimPipe produces 4 main output files:

* :ref:`First mapping BAM <output-bam>`.
* :ref:`Second mapping MAP <output-map>`.
* :ref:`Final chimeric junctions <output-chimeras>`.
* :ref:`Discarded chimeric junctions <output-chimeras>`.

.. tip:: If you want to keep intermediate output files, run ChimPipe with the ``--no-cleanup`` option. 

First mapping BAM file
~~~~~~~~~~~~~~~~~~~~~~

`BAM`_ file (``${outDir}/MappingPhase/FirstMapping/${sample_id}_firstMap.bam``) containing the reads mapped in the genome, transcriptome and *de novo* transcriptome with the `GEMtools RNA-seq pipeline`_. 

BAM is the standandard format for aligned RNA-seq reads, meaning that most analysis tools work with this format. The bam file produced can therefore be used to do other downstream analyses such as gene and transcript expression quantification.

.. _BAM: http://samtools.github.io/hts-specs/SAMv1.pdf
.. _GEMtools RNA-seq pipeline: http://gemtools.github.io/

Second mapping MAP file
~~~~~~~~~~~~~~~~~~~~~~~
MAP file (``${outDir}/MappingPhase/SecondMapping/${sample_id}_secondMap.map``) containing the reads split-mapped in the genome allowing for interchromosomal, different strand and unexpected genomic order mappings. 

Final and filtered chimeric junction files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Two tabular text files with the detected (``${outDir}/chimericJunctions_${sample_id}.txt``) and filtered out (``${outDir}/chimericJunctions_filtered_${sample_id}.txt``) chimeric splice junctions from your RNA-seq dataset. They consist on rows of 35 fields, where each row corresponds to a chimeric junction and each field contains a piece of information about the chimera. Here is a brief description of the 35 fields (most relevant fields highlighted in bold):

1. **juncCoord** - Position of the chimeric splice junction in the genome described as follows: chrA"_"coordA"_"strandA":"chrB"_"coordB"_"strandB. E. g., "chr4_90653092_+\:chr17_22023757_-" is a chimeric junction between the position 90653092 of chromosome 4 in plus strand, and the position 22023757 of chromosome chr17 in minus strand. Junction coordinates defined using 1-based system.

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


