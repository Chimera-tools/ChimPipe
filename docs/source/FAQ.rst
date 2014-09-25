.. _FAQ.rst:

====
FAQ 
====

How can I export the path to the dependencies?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To export the path of bedtools, samtools and blast (if needed) binaries you just need to type:

.. code-block:: bash

	
	$ export PATH=<BEDTOOLS_BINARIES_PATH>:<SAMTOOLS_BINARIES_PATH><BLAST_BINARIES_PATH>:$PATH
	$ # E.g. export bedtools and samtools on my system
	$ export PATH=/users/rg/brodriguez/bin/bedtools2-2.20.1/bin:/users/rg/brodriguez/bin/samtools-0.1.19:$PATH
	
How can I know the quality offset of my RNA-seq reads?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We provide a bash script to detect the offset of your RNA-seq data. You can find it at ``ChimPipe/tools/detect.fq.qual.sh``.

.. code-block:: bash

	$ bash detect.fq.qual.sh reads_1.fastq
	$ Offset 33

	
How can I know the RNA-seq library type?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We supply a bash script to infer the protocol used to produce the sequencing library from a subset of mapped read in BAM format. It is at ``ChimPipe/tools/infer_exp.sh``. You can proceed as follow:

.. code-block:: bash

	$ # 1) First you need to extract and map a subset of read pairs (E.g 1000 pairs)
	$ head -4000 reads_1.fastq > reads.subset_1.fastq
	$ head -4000 reads_2.fastq > reads.subset_2.fastq
	$ # 2) Map the reads with the gemtools rna-pipeline (at `ChimPipe/bin/gemtools-1.7.1-i3/gemtools`)	
	$ gemtools rna-pipeline -f reads.subset_1.fastq reads.subset_2.fastq -i genome.gem -a annotation.gtf -q 33 
	$ # Note: you also need your transcriptome index in the same directory as the annotation  
	$ # 3) Run our script with the BAM file generated and the annotation used for the mapping step.
	$ bash infer_exp.sh annotation.gtf reads_subset.bam
	$ bash infer_exp.sh annotation.gtf reads_subset.bam

This script will produce an output like this:

.. code-block:: bash

	Creating the reference gene model in bed format

	Inferring experiment
	Reading reference gene model gen10.long.bed ... Done
	Loading SAM/BAM file ...  Total 200000 usable reads were sampled

	This is PairEnd Data
	Fraction of reads explained by "1++,1--,2+-,2-+": 0.0135 
	Fraction of reads explained by "1+-,1-+,2++,2--": 0.9866  
	Fraction of reads explained by other combinations: 0.0000
	Removing temporary files
	DONE

The three rows starting with "Fraction of reads explained by.." are the ones giving information regarding the library. They contain several strings of three characters, i.e. 1+-, where:

* Character 1. 1 and 2 are mate1 and mate2 respectively. 
* Character 2. + and - is the strand where the read maps. 
* Character 3. + and - is the strand where the gene in which the read overlaps is annotated. 

You can apply the following rules to infer the used library from this information:

* **NONE**. Not strand-specific protocol (unstranded data). Fraction of reads explained by "1++,1--,2+-,2-+" and "1+-,1-+,2++,2--" close to 0.5000 in both cases. 

Strand-specific protocols (stranded data):
 
* **SENSE**... I am not totally sure if this protocol trully exists.
* **ANTISENSE**...  I am not totally sure if this protocol trully exists.
* **MATE1_SENSE**. Fraction of reads explained by "1++,1--,2+-,2-+" close to 1.0000. 
* **MATE2_SENSE**. Fraction of reads explained by "1+-,1-+,2++,2--" close to 1.0000.

.. note:: This script makes use of a tool from the `RSeQC package`_. So, you should have RSeQC installed in your system to run the script. We plan to recode it soon to use BEDtools instead of RSeQC and to avoid having an additional dependency.  

.. _RSeQC package: http://rseqc.sourceforge.net/.

.. tip:: In case you have any problem to map the subset of reads you can check the `GEMtools RNA Pipeline Quickstart`_ for more details. 

.. _GEMtools RNA Pipeline Quickstart: http://gemtools.github.io/docs/rna_pipeline.html


