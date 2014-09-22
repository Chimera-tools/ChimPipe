.. _FAQ.rst:

====
FAQ 
====

How can I export the path to the dependencies?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
	
	
How can I know the quality offset of my RNA-seq reads?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
	
How can I know the RNA-seq library type?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

