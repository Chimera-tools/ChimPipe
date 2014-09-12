.. _installation:

============
Installation
============

System requirements 
===================

Hardware:

* 64-bits CPU 
* ~30G RAM to run with raw RNA-Seq data from Homo sapiens

Software:

* 64-bit Linux System
* `Bedtools`_ v2.20.1 or higher  
* `Samtools`_ v0.1.19 or higher
* `Blast`_ v2.2.29+ or higher

.. _Bedtools: http://bedtools.readthedocs.org/en/latest/
.. _Samtools: http://www.htslib.org/
.. _Blast: http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download


Download ChimPipe
=================

You can do it from the `ChimPipe GitHub repository`_ in two different ways: 

.. _ChimPipe GitHub repository: https://github.com/CRG-Barcelona/ChimPipe

1. Press the `Download ZIP` button which is located at the bottom right to download the most recent version of the code as a zip archive. 
2. Clone the Git repository in case you want the code under the Git version control system (It is not allowed to contribute to the project, but it will be in the future):

.. code-block:: bash

	# Move into the folder you want to clone the repositoy.
	$ cd /users/rg/brodriguez/bin
	# Clone it.
	$ git clone git@github.com:brguez/ChimPipe.git

Once you downloaded it, if you go to the directory you will see that it is organized in the following way:

.. code-block:: bash


Set up ChimPipe's path
======================
Finally, go to your newly created ChimPipe's directory and do make to set the path to ChimPipe in your system. 

.. code-block:: bash

	# Move into ChimPipe.
	$ cd /users/rg/brodriguez/bin/ChimPipe
	# Do make to set the path to ChimPipe in the system 
	$ make
	
You should get the following message if everything goes well:

.. code-block:: bash

	1) Make a temporary bash script to write in ChimPipe.sh the code to define a variable 
	containing the path to the directory where the pipeline is placed... done
	2) Execute the script... done 
	3) Remove the script... done
	

