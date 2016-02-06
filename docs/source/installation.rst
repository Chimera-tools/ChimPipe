.. _installation:

============
Installation
============

Download
=========

You can do it from our `GitHub repository`_ in two different ways: 

.. _GitHub repository: https://github.com/Chimera-tools/ChimPipe.git

1. Go to the releases tab and download the latest release.    
2. Clone the git repository in case you want the latest version of the code:

.. code-block:: bash

	# Move into the folder in which you want to clone the repositoy.
	$ cd ~/apps
	# Clone it.
	$ git clone https://github.com/Chimera-tools/ChimPipe.git

ChimPipe does not require any further installation step. It already comes with precompiled gemtools binaries. It is written in Bash and Awk and can be run as a standalone application on diverse Linux systems.

Requirements 
=============

Hardware:

* 64-bits CPU. 
* RAM: ~40G for 100 million illumina paired-end reads.

Software:

* 64-bit Linux System.
* `Bedtools`_ v2.20.1 or higher.  
* `Samtools`_ v0.1.19 or higher.
* `Blast`_ v2.2.29+ or higher (only if you want to generate your own similarity text files, see **Gene pair similarity file** in the :ref:`manual` section).

.. _Bedtools: http://bedtools.readthedocs.org/en/latest/
.. _Samtools: http://www.htslib.org/
.. _Blast: http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download








