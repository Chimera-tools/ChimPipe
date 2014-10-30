.. _installation:

============
Installation
============

System's requirements 
======================

Hardware:

* 64-bits CPU 
* RAM: ~40G for 100 million human PE reads of 52 bases. 

Software:

* 64-bit Linux System
* `Bedtools`_ v2.20.1 or higher  
* `Samtools`_ v0.1.19 or higher
* `Blast`_ v2.2.29+ or higher (only if you want to generate your own similarity text files, see **Reference between gene pairs** in the :ref:`manual` section)

.. _Bedtools: http://bedtools.readthedocs.org/en/latest/
.. _Samtools: http://www.htslib.org/
.. _Blast: http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download


Downloading ChimPipe
====================

You can do it from the `ChimPipe GitHub repository`_ in two different ways: 

.. _ChimPipe GitHub repository: https://github.com/Chimera-tools/ChimPipe.git

1. Press the `Download ZIP` button which is located at the bottom right of the page to download the most recent version of the code as a zip archive. 
2. Clone the Git repository in case you want the code under the Git version control system (currently, it is not allowed to contribute to the project, but it will in the future):

.. code-block:: bash

	# Move into the folder in which you want to clone the repositoy.
	$ cd /users/rg/brodriguez/bin
	# Clone it.
	$ git clone https://github.com/Chimera-tools/ChimPipe.git


Setting up ChimPipe's path
==========================
Finally, go to your newly created ChimPipe's directory and type make to set up the path to ChimPipe in your system. 

.. code-block:: bash

	# Move into ChimPipe.
	$ cd /users/rg/brodriguez/bin/ChimPipe
	# Type make to set up the path to ChimPipe in your system 
	$ make
	
	

