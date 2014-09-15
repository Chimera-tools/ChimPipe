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
	

* Illumina CASAVA package 1.8 and higher. The identifier consists on two strings separated by a space with the second one specifying the mate in the first digit. E. g:   

.. code-block:: bash
	
	$ # Mate 1
	@seq.1_set1:140:D034KACXX:1:2105:18345:62400 1:N:0:
	CCCAGCCTGTTTCCCTGCCTGGGAAACTAGAAAGAGGGCTTCTTTCTCCT
	+
	IJJJIIIIIGIGIEBHHHGFFFF=CDEEEEEDDDDCCCDDA?55><CBBB
	
	$ # Mate 2
	@seq.1_set1:140:D034KACXX:1:2105:18345:62400 2:N:0:
	GCACCCTTCACTCCCTCCCTTGGGCGCCTCCCTCCCGAGGGTAGGGACCC
	+
	FFHIJJCHIIHBHIIIAHFFFFFCDEDEECDBB;??@CD?CCCCCCC@CC

I case your FASTQ files do not meet this format you should modify the identifier. Awk is a perfect tool for such kind of problems. E. g: I downloaded a dataset from the NCBI Sequence Read archive:

.. code-block:: bash
	
	$ cd /users/rg/brodriguez/bin/ChimPipe/reads/berger
	
	$ # But it does not have a proper identifier
	$ head -4 berger_1.fastq
	
	@SRR018259.1 BI:080831_SL-XAN_0004_30BV1AAXX:5:1:708:1868 length=51
	GTAACATATTCACAGACATGTCAGTGTGCACTCAGGAACACTGGTTTCATT
	+SRR018259.1 BI:080831_SL-XAN_0004_30BV1AAXX:5:1:708:1868 length=51
	IIIIIIIIIIIIIII>=CII=8=H032/-++D+'.@)2+4/+)1'4.#"*.
	
	$ ## So, I use awk to modify it. 
	
	$ awk '{if (NR%2==1){print $1"_"$2"_"$3"#0/1"} else {print $0}} berger_1.fastq > berger_fixed_1.fastq		
	$ head -4 berger_fixed_1.fastq 
	@SRR018259.1_BI:080831_SL-XAN_0004_30BV1AAXX:5:1:708:1868_length=51#0/1
	GTAACATATTCACAGACATGTCAGTGTGCACTCAGGAACACTGGTTTCATT
	+SRR018259.1_BI:080831_SL-XAN_0004_30BV1AAXX:5:1:708:1868_length=51#0/1
	IIIIIIIIIIIIIII>=CII=8=H032/-++D+'.@)2+4/+)1'4.#"*.

	$ # Finally, I apply the same procedure for the mate 2..

Genome index
~~~~~~~~~~~~
An indexed reference genome in GEM format has to be provided to do the mapping steps. You just need to run the GEMtools indexer (supplied with ChimPipe distribution) with your genome in FASTA format to produce it:

.. _FASTA:
 
.. code-block:: bash

	$ cd /users/rg/brodriguez/bin/ChimPipe/genomes/GRCh37
	$ gemtools=/users/rg/brodriguez/bin/ChimPipe/bin/gemtools-1.7.1-i3/gemtools
	$ $gemtools index -i Homo_sapiens.GRCh37.chromosomes.chr.M.fa

Note that you can specify multiple threads with the option -t. You should get the following message if everything goes well:

.. code-block:: bash


Genome annotation
~~~~~~~~~~~~~~~~~ 
Chimpipe also takes as input a genome annotation in `GTF`_ format to find reads spanning splice junctions between exons from two different genes. This annotation has to contain at least one tag-value pair in the attributes field with the gene id and two optional pairs will be taken into account by ChimPipe if supplied: gene name and gene type. E.g:

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
As explained in the :ref:`installation` section, you need to have installed BEDtools and SAMtools to execute ChimPipe, plus blast in case you want to produce your own similarity between gene pairs text files (See **Similarity between gene pairs**). In case you do not have them, you can not download an install them from their webpages. Once installed, you have to export the path to their binaries as follow:  

.. code-block:: bash

	# Check the current path setting
	$ echo "$PATH"
	$ /software/rg/el6.3/pythonz/bin:/software/rg/el6.3/bin:/software/rg/el6.3/texlive/2012/bin/x86_64-linux:/software/as/el6.3/test/modules/Modules/3.2.10/bin/:/usr/lib64/qt-3.3/bin:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/usr/lib64/openmpi/bin/:/usr/lib64/compat-openmpi/bin/:/users/rg/brodriguez/bin:/software/rg/bin/

	# Add the path to bedtools, samtools and blast binaries to the $PATH variable
	$ export PATH=/users/rg/brodriguez/bin/bedtools2-2.20.1/bin:/users/rg/brodriguez/bin/samtools-0.1.19:/users/rg/brodriguez/bin/blast-2.2.29+/bin:$PATH
	$ echo $PATH
	$ /users/rg/brodriguez/bin/bedtools2-2.20.1/bin:/users/rg/brodriguez/bin/samtools-0.1.19:/users/rg/brodriguez/bin/blast-2.2.29+/bin:/software/rg/el6.3/pythonz/bin:/software/rg/el6.3/bin:/software/rg/el6.3/texlive/2012/bin/x86_64-linux:/software/as/el6.3/test/modules/Modules/3.2.10/bin/:/usr/lib64/qt-3.3/bin:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/usr/lib64/openmpi/bin/:/usr/lib64/compat-openmpi/bin/:/users/rg/brodriguez/bin:/software/rg/bin/
	# Now, if you can call directly bedtools, samtools or blast to check if it is working. E.g:
	$ samtools
	

2. Check the quality offset (Skip if you already know)  
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The quality scores (Q) measure the probability that a base is called incorrectly by the sequencing machine. Within your FASTQ files, they are represented in the fourth line of each read as an string of ASCII characters (each character correspond to the Q score of a certain base in the sequencing read). The correspondence between each ASCII character and the Q score is based on some offset. These offset vary with the sequencing platform (current Illumina machines uses 33, while older ones 33). 

Along with ChimPipe, we supply an in-house bash script you can use to detect the offset of your reads before running ChimPipe:

.. code-block:: bash

	$ cd /users/rg/brodriguez/bin/ChimPipe/reads/berger
	$ ls 
	$ berger_1.fastq berger_2.fastq
	$ # I will use the detect.fq.qual to know the offset. 
	$ quality=/users/rg/brodriguez/bin/ChimPipe/tools/detect.fq.qual.sh
	$ $quality berger_1.fastq
	$ Offset 33
	$ # Ok, the offset is 33. I will use this information to run the pipeline afterwards. 

3. Check the RNA-seq library type (Skip if you already know)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Different protocols that can be used to generate a RNA-seq library. There are also important differences among them that have to be taken into account in several steps of the chimera detection pipeline. However, ChimPipe can not determine the protocol used to produce your reads, so you need to supply this information manually with the option **--read-directionality <STRING>**. Where **STRING** has to be one of these library types:

* NONE. Not strand-specific protocol (unstranded data). The information about from which strand the transcript is transcribed is not available. **Default configuration**

Strand-specific protocols (stranded data):
 
* SENSE. Transcript directly sequenced. Reads map to the transcript strand.
* ANTISENSE. Reverse complementary sequence of the transcript sequenced. Reads map to the opposite strand of the transcript. 
* MATE1_SENSE. Reads on the left of the fragment (mates 1) sequenced from the transcript (map in the transcrip strand), and the ones in the right (mates 2) sequenced from the complementary reverse sequence (map in the opposite strand). 
* MATE2_SENSE. Reads on the left of the fragment (mates 1) sequenced from the complementary reverse sequence (map in the opposite strand), and the ones in the right (mates 2) sequenced from the transcript (map in the transcrip strand). 

In case you do not know the protocol used to produce your data, you can use a tool provided along ChimPipe to infer it from a subset of mapped reads as follows:

.. code-block:: bash

	$ cd /users/rg/brodriguez/bin/ChimPipe/reads/berger
	$ # I will substract a subset of 1000 reads 
	$ head -4000 berger_1.fastq > berger.subset1000_1.fastq
	$ head -4000 berger_2.fastq > berger.subset1000_2.fastq
	$ gemtools=/users/rg/brodriguez/bin/ChimPipe/bin/gemtools-1.7.1-i3/gemtools
	$ index=
	$ annot=
	$ gemtools --loglevel $loglevel rna-pipeline -f berger.subset1000_1.fastq -i $index -a $annot -q  -n berger 
	

4. Run ChimPipe
~~~~~~~~~~~~~~~

Output
======

By default, ChimPipe produces 3 files as output:

* First mapping BAM file
* Second mapping MAP file
* Chimeric junctions text file



First mapping BAM file
~~~~~~~~~~~~~~~~~~~~~~
`BAM`_ file containing the reads mapped in the genome, transcriptome and *de novo* transcriptome with the GEMtools RNA pipeline. 

Many next-generation sequencing analysis tools work with this format, so it can be used to do very different analyses such as gene and transcript quantification or differential gene expression analysis.

.. _BAM: http://samtools.github.io/hts-specs/SAMv1.pdf

Second mapping MAP file
~~~~~~~~~~~~~~~~~~~~~~~
MAP file containing reads segmentally mapped in the genome allowing for interchromosomal, different strand and unexpected genomic order mappings. 

*Unfortunatelly, there is not an official documentation describing this mapping format. 

Chimeric junctions text file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
chimeric_junctions.txt

Here is a brief description of the GFF fields:

1. **juncId** -  	
2. **nbstag** - 
3. **nbtotal** -
4. **maxbeg** -
5. **maxEnd** - 
6. **samechr** -
7. **samestr** -
8. **dist** -
9. **ss1** -
10. **ss2**	- 
11. **gnlist1** -	
12. **gnlist2**	- 
13. **gnname1** - 
14. **gnname2**	- 
15. **bt1** -	
16. **bt2**	- 
17. **PEsupport** -
18. **maxSim** -	
19. **maxLgal** -

**Example**

.. code-block:: bash

	Here's an example of a GFF-based track.

	juncId	nbstag	nbtotal	maxbeg	maxEnd	samechr	samestr	dist	ss1	ss2	gnlist1	gnlist2	gnname1	gnname2	bt1	bt2	PEsupport	maxSim	maxLgal
	chr1_121115975_+:chr1_206566046_+ 1 1 121115953 206566073 1 1 85450071 GC AG SRGAP2D, SRGAP2,SRGAP2C, SRGAP2D, SRGAP2,SRGAP2C, . . 1-1:2,1-2:2, 99.44 1067

