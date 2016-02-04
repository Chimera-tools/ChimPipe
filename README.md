## Welcome to ChimPipe's repository!

ChimPipe is a pipeline for the detection of chimeric transcripts and gene fusions from Paired-End RNA-seq data developed at the Centre for Genomic Regulation (CRG) in Barcelona. 

#### Download ChimPipe
Two different ways:

* Go to the releases tab and download the latest release. 

* Clone ChimPipe's Git repository in case you want the latest version of the code (currently, it is not allowed to contribute to the project):

```
# Move into the folder in which you want to clone the repositoy.
$ cd ~/apps
# Clone it.
$ git clone https://github.com/Chimera-tools/ChimPipe.git
```

ChimPipe does not require any installation step. It already comes with precompiled gemtools binaries. It is written in Bash and Awk and can be run as a standalone application on diverse Linux systems. 

#### System’s requirements

1. Hardware:

    * 64-bits CPU
    * RAM: ~40G for 100 million illumina paired-end reads.

2. Software:

    * 64-bit Linux System
    * Bedtools v2.20.1 or higher
    * Samtools v0.1.19 or higher
    * Blast v2.2.29+ or higher 

#### ChimPipe's documentation
Please check [ChimPipe's documentation] (https://chimpipe.readthedocs.org/) to find detailed information about:

* How ChimPipe works
* Install the pipeline
* Run the pipeline
* Contact us


