

How can I compute the gene pair similarity file?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
One of ChimPipe's steps to discard actefactual chimeras is to filter out those chimeric junctions connecting genes that encode transcripts with a high sequence homology. Although it is an optional filter, it is **strongly recommended**, since according to our benchmark it improves a lot the specificity with a minimal decrease on the sensitivity.  

To enable this homology-based filtering, you only need to run ChimPipe with the option ``--similarity-gene-pairs <TEXT FILE>``, where **TEXT FILE** is a file containing the matrix with information about the sequence similarity between gene pairs. 

You can download our **pre-generated matrices** por human, mouse and drosophila annotations from :ref:`Downloads` section or you can produce your own matrix with the script ``ChimPipe/src/bash/tools/similarity_bt_gnpairs.sh`` as follows:

	$ bash similarity_bt_gnpairs.sh annot.gtf genome.gem

Please check out our :ref:`FAQ <faq-similarity>` section for more information about how the script works. 
	
.. warning:: Make sure you run ChimPipe with a similarity matrix generated from the annotation and genome you are using.  



.. tip:: It is strongly advisable to use the option ``--similarity-gene-pairs <TEXT FILE>`` to discard junctions connecting genes which encode transcripts with a high sequence homology, which are likely sequencing or mapping artefacts. Please check :ref:`Gene pair similarity file <input-similarity>` section above to learn how to produce the text file needed. 

