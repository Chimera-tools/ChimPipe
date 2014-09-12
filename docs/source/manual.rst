.. _manual:

======
Manual
======

This section describes the files ChimPipe takes as input, how to generate them, run the pipeline and interpret its output. 

Input files
===========
ChimPipe needs 4 mandatory different types of files plus 1 optional.  

Mandatory:

* Paired-end (PE) RNA-seq reads
* Genome index 
* Transcriptome annotation
* Transcriptome annotation index

Optional:

* Similarity between gene pairs text file


Paired-end (PE) RNA-seq reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ChimPipe has been designed to deal with `Illumina paired-end`_ RNA sequencing data. It takes two `FASTQ`_ files as input, one for the first mates and another one for the second mates in the read pairs respectively. Make sure that they are in the same directory and they are named according to this convention: 

.. _Illumina paired-end: http://technology.illumina.com/technology/next-generation-sequencing/paired-end-sequencing_assay.ilmn
.. _FASTQ: http://maq.sourceforge.net/fastq.shtml

* **Mate 1**. "SampleId + [.1|_1] + [.fastq|.txt] (+ .gz if compressed)"
* **Mate 2**. "SampleId + [.2|_2] + [.fastq|.txt] (+ .gz if compressed)"

E. g. BERGER_1.fastq.gz and BERGER_2.fastq.gz would be the compressed FASTQ files for mate 1 and mate 2 in the BERGER sample. Note that you should use the same convention for both mates, so if mate 1 is BERGER_1.fastq.gz, mate 2 can not be BERGER.2.txt.gz

The FASTQ file uses four lines per sequencing read. You need to check the format of the first line of each read, which begins with the '@' character and is followed by an identifier. This identifier should meet one of the two Illumina standards to specify which member of the pair the read is:

* Illumina CASAVA package lower than 1.8. The identifier has to be a single string ended in /1 or /2 for mate 1 and mate 2, respectively. E. g.:

.. code-block:: bash
	
	# Mate 1
	@SRR018259.1_BI:080831_SL-XAN_0004_30BV1AAXX:5:1:708:1868#0/1
	GTAACATATTCACAGACATGTCAGTGTGCACTCAGGAACACTGGTTTCATT
	+SRR018259.1_BI:080831_SL-XAN_0004_30BV1AAXX:5:1:708:1868#0/1
	IIIIIIIIIIIIIII>=CII=8=H032/-++D+'.@)2+4/+)1'4.#"*.
	
	# Mate 2
	@SRR018259.1_BI:080831_SL-XAN_0004_30BV1AAXX:5:1:708:1868#0/2
	CAGATGGCATTGCAGAACTAAAAAAGCAAGCTGAATCAGTGTTAGATCTCC
	+SRR018259.1_BI:080831_SL-XAN_0004_30BV1AAXX:5:1:708:1868#0/2
	IIIGIIIIIIIIIFI:>DID<AFH=>78@I41I055549.?+.42-'%**'
	

* Illumina CASAVA package 1.8 and higher. The identifier consists on two strings separated by a space with the second one specifying the mate in the first digit. E. g:   

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
	
	$ cd /users/rg/brodriguez/bin/ChimPipe/reads/berger
	
	# But it does not have a proper identifier
	$ head -4 berger_1.fastq
	
	@SRR018259.1 BI:080831_SL-XAN_0004_30BV1AAXX:5:1:708:1868 length=51
	GTAACATATTCACAGACATGTCAGTGTGCACTCAGGAACACTGGTTTCATT
	+SRR018259.1 BI:080831_SL-XAN_0004_30BV1AAXX:5:1:708:1868 length=51
	IIIIIIIIIIIIIII>=CII=8=H032/-++D+'.@)2+4/+)1'4.#"*.
	
	## So, I use awk to modify it. 
	
	$ awk '{if (NR%2==1){print $1"_"$2"_"$3"#0/1"} else {print $0}} berger_1.fastq > berger_fixed_1.fastq		
	$ head -4 berger_fixed_1.fastq 
	@SRR018259.1_BI:080831_SL-XAN_0004_30BV1AAXX:5:1:708:1868_length=51#0/1
	GTAACATATTCACAGACATGTCAGTGTGCACTCAGGAACACTGGTTTCATT
	+SRR018259.1_BI:080831_SL-XAN_0004_30BV1AAXX:5:1:708:1868_length=51#0/1
	IIIIIIIIIIIIIII>=CII=8=H032/-++D+'.@)2+4/+)1'4.#"*.

	# Finally, I apply the same procedure for the mate 2..

Genome index
~~~~~~~~~~~~
An indexed reference genome in GEM format has to be provided to do the mapping steps. You just need to run the GEMtools indexer (provided with ChimPipe distribution) with your genome in FASTA format to produce it:

.. _FASTA:
 
.. code-block:: bash

	$ cd /users/rg/brodriguez/bin/ChimPipe/genomes/GRCh37
	$ gemtools=/users/rg/brodriguez/bin/ChimPipe/bin/gemtools-1.7.1-i3/gemtools
	$ $gemtools index -i Homo_sapiens.GRCh37.chromosomes.chr.M.fa

Note that you can specify multiple threads with the option -t. You should get the following message if everything goes well:

.. code-block:: bash


Genome annotation
~~~~~~~~~~~~~~~~~ 
Chimpipe also takes as input a genome annotation in `GTF`_ format to find reads spanning splice junctions between exons from two different genes. This annotation has to contain at least one tag-value pair in the attributes field with the gene id and two optional pairs will be taken into account by ChimPipe if provided: gene name and gene type. E.g:

.. _GTF: http://www.ensembl.org/info/website/upload/gff.html



.. code-block:: bash
	
	# This is an example of one annotated exon with an appropiated format. 	
	# The attributes are the gene id (mandatory), the gene type and gene name (optional), 
	# plus some additional tag-value pairs that will not be considered by ChimPipe.   
	
	chr1	HAVANA	exon	69091	70008	.	+	.	gene_id "ENSG00000186092.4"; transcript_id "ENST00000335137.3"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "OR4F5";
	transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "OR4F5-001"; exon_number 1; exon_id "ENSE00002319515.1"; level 2; tag "basic"; tag "appris_principal";
	tag	"CCDS"; ccdsid "CCDS30547.1"; havana_gene "OTTHUMG00000001094.1"; havana_transcript "OTTHUMT00000003223.1";

ChimPipe has been benchmarked with `Gencode v10`_ and `UCSC Known Genes`_ annotation. It displayed a better sensitivity with Gencode v10 while the similar false positive rate was similar (see Benchmark section). Thus, we encourage the user to use Gencode annotation, it is a richer annotation what increase the sensitivity of the chimera detection process. 

.. _Gencode v10: http://www.gencodegenes.org/releases/10.html
.. _UCSC Known Genes: https://genome.ucsc.edu/cgi-bin/hgTables?command=start

Transcriptome annotation index
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
An indexed transcriptome annotation in GEM format has to be given as input to find reads spanning annotated splice junctions. You only have to run the GEMtools transcriptome indexer (provided with ChimPipe distribution) with your GEM indexed genome and its annotation in GTF format to generate it. 

.. code-block:: bash

	$ cd /users/rg/brodriguez/bin/ChimPipe/annotations/gencode10
	$ gemtools=/users/rg/brodriguez/bin/ChimPipe/bin/gemtools-1.7.1-i3/gemtools
	$ genome=/users/rg/brodriguez/bin/ChimPipe/genomes/GRCh37/Homo_sapiens.GRCh37.chromosomes.chr.gem
	$ $gemtools t-index -i $genome -a gen10.long.gtf	

You can specify multiple threads with -t. You should get the following message if everything goes well: 

.. code-block:: bash

**IMPORTANT**: The indexed gene annotation has to be placed in the same folder as the genome annotation to be used by ChimPipe

Similarity between gene pairs (Optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Execute pipeline
================

1. Set up the environment
~~~~~~~~~~~~~~~~~~~~~~~~~
As explained in the installation section, Export the paths of bedtools and samtools

2. Check the quality offset of your reads  
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Illumina uses Phred quality scores. A Phred quality score is produced by a model that uses quality predictors as inputs and produces Q scores as outputs. Illumina's algorithm uses many predictors to generate a Q score for each base in every read.. These predictors can be divided into three categories: single base, full-read and “neighborhood”.

A Q score is simply a representation of the probability that a base-call is incorrect. Here's the equation: Q = -10 * log10p where p is the probability that the particular base-call is incorrect (calculated using the set of predictors). A base-call with a Q score of 20 has a 1% probability of being incorrect; Q30 corresponds to 0.1%, Q40 to 0.01%, and so on.

Note about offsets:

Within your sequencing files Q scores are generally represented using a single ASCII character. Each ASCII character has a value that corresponds to the Q score value based on some offset. Sanger Q scores use 33 as the offset value. Currently Illumina also uses 33, but that wasn't always the case. Older Illumina data used 64 as an offset. This isn't really a big deal. The quality of the data isn't affected, but some secondary analysis tools will need to be told which offset has been used. You can find this information in the header of the Per Base Sequence Quality graph. 

.. code-block:: bash



3. Check the directionality of your reads 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

4. Run ChimPipe
~~~~~~~~~~~~~~~

Output
======

ChimPipe produces 3 files as output:

* First mapping BAM file
* Second mapping MAP file
* Chimeric junctions text file

First mapping BAM file
~~~~~~~~~~~~~~~~~~~~~~


Second mapping MAP file
~~~~~~~~~~~~~~~~~~~~~~~


Chimeric junctions text file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


