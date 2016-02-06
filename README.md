## ChimPipe

ChimPipe is a computational method for the detection of novel chimeric transcripts and fusion genes from Illumina Paired-End RNA-seq data. It combines junction spanning and paired-end read information to accurately detect chimeric splice junctions at base-pair resolution. ChimPipe have been developed at the Centre for Genomic Regulation (CRG) in Barcelona. 

## Download 
Two different ways:

* Go to the releases tab and download the latest release. 

* Clone ChimPipe's Git repository in case you want the latest version of the code:

```
# Move into the folder in which you want to clone the repositoy.
$ cd ~/apps
# Clone it.
$ git clone https://github.com/Chimera-tools/ChimPipe.git
```

ChimPipe does not require any further installation step. It already comes with precompiled GEMtools binaries. It is written in Bash and Awk and can be run as a standalone application on diverse Linux systems. 

## Requirements

1. Hardware:

    * 64-bits CPU
    * RAM: ~40G for 100 million illumina paired-end reads.

2. Software:

    * 64-bit Linux System
    * Bedtools v2.20.1 or higher
    * Samtools v0.1.19 or higher
    * Blast v2.2.29+ or higher 

## Documentation
Please check [ChimPipe's documentation] (https://chimpipe.readthedocs.org/) to find detailed information about:

* How ChimPipe works
* Install ChimPipe
* Execute ChimPipe
* Contact with us

## License
ChimPipe is distributed under the GPLv3. Consult the [LICENSE] (https://github.com/Chimera-tools/ChimPipe/blob/master/LICENSE.txt) file for more information.


