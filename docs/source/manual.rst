.. _manual:

======
Manual
======

This section describes ChimPipe' input files, how to generate them, how to run ChimPipe, and how to interpret its outputs. 

Input files
===========
ChimPipe takes as input 4 mandatory files and 1 optional file.  

Mandatory:

* :ref:`Paired-end (PE) RNA-seq reads <input-reads>`
* :ref:`Genome index <input-genome-index>` 
* :ref:`Genome annotation <input-annot>`
* :ref:`Transcriptome index <input-annot-index>`

Optional:

* :ref:`Gene pair similarity file <input-similarity>`

.. _input-reads:

Paired-end (PE) RNA-seq reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ChimPipe has been designed to deal with `Illumina paired-end`_ RNA sequencing data. It takes as input two `FASTQ`_ files, one for the first mates and one for the second mates of the read pairs. These files need to be located in the same directory and named according to this convention: 

.. _Illumina paired-end: http://technology.illumina.com/technology/next-generation-sequencing/paired-end-sequencing_assay.ilmn
.. _FASTQ: http://maq.sourceforge.net/fastq.shtml

* **Mate 1**. "SampleId + [.1|_1] + [.fastq|.txt] (+ .gz if compressed)"
* **Mate 2**. "SampleId + [.2|_2] + [.fastq|.txt] (+ .gz if compressed)"

E. g. BERGER_1.fastq.gz and BERGER_2.fastq.gz would be the compressed FASTQ files for mate 1 and mate 2 in the BERGER sample. 

.. warning:: The **read identifier** should follow **Illumina convention** to specify which member of the pair the read is. See our :ref:`FAQ <faq-reads>` section for more information. 

.. warning:: Make sure the 2 FASTQ files are in the **same directory** and you use the **same convention** for both mates. E. g: if mate 1 is BERGER_1.fastq.gz, mate 2 can not be BERGER.2.txt.gz but has to be BERGER_2.fastq.gz

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

Genome annotation
~~~~~~~~~~~~~~~~~
ChimPipe also takes as input a genome annotation file in `GTF`_ format (including annotated exons). It can contain other features different from exons, i. e. introns or UTR, but they will not be considered by ChimPipe in the chimera detection process. This annotation needs to contain at least two (tag,value) pairs in the attribute field with the gene_id and the transcript_id tags. Two optional (tag,value) pairs will be taken into account by ChimPipe if they are provided: gene_name and gene_type. e.g:

.. _GTF: http://www.ensembl.org/info/website/upload/gff.html

.. code-block:: bash
	
	# This is an example of an annotated exon with an appropiate format. 	
	# The attributes are the gene_id, transcript_id (mandatory), the gene type and gene name (optional), 
	# plus some additional (tag,value) pairs that will not be considered by ChimPipe.   
	
	chr1	HAVANA	exon	69091	70008	.	+	.	gene_id "ENSG00000186092.4"; transcript_id "ENST00000335137.3"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "OR4F5";
	transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "OR4F5-001"; exon_number 1; exon_id "ENSE00002319515.1"; level 2; tag "basic"; tag "appris_principal";
	tag	"CCDS"; ccdsid "CCDS30547.1"; havana_gene "OTTHUMG00000001094.1"; havana_transcript "OTTHUMT00000003223.1";

.. note:: ChimPipe has been benchmarked with `Gencode v10`_ and `UCSC Known Genes`_  human gene annotations. It displayed a better sensitivity with Gencode v10 but a similar false positive rate with both annotations. Thus, it is advisable to use Gencode annotation, since it is a richer annotation which increases the sensitivity of the chimera detection process. 

.. _Gencode v10: http://www.gencodegenes.org/releases/10.html
.. _UCSC Known Genes: https://genome.ucsc.edu/cgi-bin/hgTables?command=start

.. _input-annot-index:

Transcriptome index
~~~~~~~~~~~~~~~~~~~
An transcriptome index in GEM format has to be provided as input to ChimPipe in the same directory as the genome annotation GTF file, in order to find reads spanning annotated splice junctions. 

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

.. warning:: The transcriptome index has to be placed in the same folder as the genome annotation to be used by ChimPipe.

.. _input-similarity:


Executing ChimPipe
==================

1. Setting up the environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
As explained in the :ref:`installation` section, you need to have BEDtools, SAMtools and Blast installed on your system to execute ChimPipe. In case you do not have them, you can download and install them from their web pages. Once they are installed, you have to export the path to their binaries. 

Please check out our :ref:`FAQ <faq-dependencies>` section in case you have any problem.  
	
2. Running ChimPipe
~~~~~~~~~~~~~~~~~~~
Once you have generated the genome and the transcriptome indices, you know the quality offset and the library type of your PE RNA-seq reads, you can run ChimPipe as follows:

.. code-block:: bash
	
	bash ChimPipe.sh -i reads_1.fastq -g genome.gem -a annotation.gtf 

All these files and parameters given as input to ChimPipe are **mandatory arguments**. Please see bellow their descripion: 

.. code-block:: bash

	-i|--input reads_1.fastq – First mate sequencing reads. ChimPipe deals with paired-end data. 
				   Please make sure the second mate file is in the same directory as 
				   the first one, and the files are named according to the same convention. 
				   E.g: the second mate of "reads_1.fastq" should be "reads_2.fastq". 
						   
	-g|--genome-index genome.gem – Index for the reference genome in GEM format.

	-a|--annotation annotation.gtf – Reference genome annotation file in GTF format. The transcriptome 
						index has to be in the same directory as the annotation. 
								 
**Optional arguments.** Please do ``ChimPipe.sh -h or --help`` to see a short help with the main options. You can also do ``ChimPipe.sh --full-help`` to see the all the possible options. 

.. tip:: If your machine has more than one CPU it is recommended to run ChimPipe with multiple threads (at least 4). It will speed up the mapping steps a lot. Use the option ``-t|--threads <threads>``, where **threads** is the number of CPUs available. 

.. note:: The pipeline is restartable: this means that if ChimPipe fails and you run it again, it will skip the steps that have been already completed. You just need to make sure that you removed the files generated in the step where the pipeline failed. 

Output
======

By default, ChimPipe produces 3 output files:

* :ref:`First mapping BAM file <output-bam>` 
* :ref:`Second mapping MAP file <output-map>` 
* :ref:`Chimeric junctions file <output-chimeras>` 

.. tip:: If you want to keep intermediate output files, run ChimPipe with the ``--no-cleanup`` option. 

.. _output-bam:

First mapping BAM file
~~~~~~~~~~~~~~~~~~~~~~
`BAM`_ file containing the reads mapped in the genome, transcriptome and *de novo* transcriptome with the `GEMtools RNA-seq pipeline`_. 

This is the standard format for next-generation sequencing, meaning that most analysis tools work with this format. The bam file produced can therefore be used to do other downstream analyses such as gene and transcript quantification or differential gene expression analysis.

.. _BAM: http://samtools.github.io/hts-specs/SAMv1.pdf
.. _GEMtools RNA-seq pipeline: http://gemtools.github.io/

.. _output-map:

Second mapping MAP file
~~~~~~~~~~~~~~~~~~~~~~~
MAP file containing reads segmentally mapped in the genome allowing for interchromosomal, different strand and unexpected genomic order mappings. 

.. _output-chimeras:

Chimeric junction file
~~~~~~~~~~~~~~~~~~~~~~
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

