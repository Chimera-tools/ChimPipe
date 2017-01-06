## ChimPipe

ChimPipe is a computational method for the detection of novel transcription-induced chimeric transcripts and fusion genes from Illumina Paired-End RNA-seq data. It combines junction spanning and paired-end read information to accurately detect chimeric splice junctions at base-pair resolution. 

ChimPipe have been developed at the [Computational Biology of RNA Processing group] (http://www.crg.eu/en/roderic_guigo) in Barcelona, Spain.

## Download 
Two different ways:

* Go to the releases tab and download the latest release. 

* Clone the git repository in case you want the latest version of the code:

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
Please check ChimPipe [documentation] (https://chimpipe.readthedocs.org/) to find detailed information about:

* How ChimPipe works
* Installation
* Execute ChimPipe

## License
ChimPipe is distributed under the GPLv3. Consult the [LICENSE] (https://github.com/Chimera-tools/ChimPipe/blob/master/LICENSE.txt) file for more information.

## Contact
Please feel free to join the ChimPipe’s [mailing list] (https://groups.google.com/forum/#!forum/chimpipe-mailing-list/) in case you have a question, you want to report an issue or request a feature.

You can also directly contact us at: chimpipe.pipeline@gmail.com


## Citation

Please cite the following article if you use ChimPipe in your research:

```
Rodríguez-Martín B, Palumbo E, Marco-Sola S, Griebel T, Ribeca P, Alonso G, Rastrojo A, Aguado B, 
Guigó R, Djebali S. (2017) ChimPipe: accurate detection of fusion genes and transcription-induced 
chimeras from RNA-seq data. BMC Genomics 18(1):7
```


