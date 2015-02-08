###Welcome to ChimPipe's repository!

ChimPipe is a pipeline for the detection of chimeric transcripts from Paired-End RNA-seq data developed at the Centre for Genomic Regulation (CRG) in Barcelona. 


## Download ChimPipe
Two different ways:

* Press the Download ZIP button which is located at the bottom right of the page to download the most recent version of the code as a zip archive.

* Clone the Git repository in case you want the code under the Git version control system (currently, it is not allowed to contribute to the project, but it will in the future):

```
# Move into the folder in which you want to clone the repositoy.
$ cd /users/rg/brodriguez/bin
# Clone it.
$ git clone https://github.com/Chimera-tools/ChimPipe.git
```

## Install ChimPipe
Go to your newly created ChimPipe’s directory and type `make` to set up the path to ChimPipe in your system.

## System’s requirements

1. Hardware:

    * 64-bits CPU
    * RAM: ~40G for 100 million human PE reads of 52 bases.

2. Software:

    * 64-bit Linux System
    * Bedtools v2.20.1 or higher
    * Samtools v0.1.19 or higher
    * Blast v2.2.29+ or higher (only if you want to generate your own similarity text files, see Reference between gene pairs in the Manual section)

## ChimPipe's documentation
Please check [ChimPipe's documentation] (https://chimpipe.readthedocs.org/) to find detailed information about:

* How ChimPipe works
* Install the pipeline
* Run the pipeline
* Contact with us


