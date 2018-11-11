# Zika virus E protein deep mutational scanning
Deep mutational scanning of Zika virus E protein.
Experiments performed by the [Matt Evans lab](http://labs.icahn.mssm.edu/evanslab/).
Sequencing and computational analyses performed by the [Bloom lab](https://research.fhcrc.org/bloom/en.html).

## Analysis and results
The analysis is performed by the Jupyter notebook [analysis_notebook.ipynb](analysis_notebook.ipynb) using [dms_tools2](https://jbloomlab.github.io/dms_tools2/).
That notebook also contains the results.

## Input data
The input data are in the [./data/](data) subdirectory. 
These data consist of:

 - [./data/E.fasta](data/E.fasta): coding sequence of E protein from ZIKV MR766 strain.

 - [./data/samplelist.csv](data/samplelist.csv): all the samples that we sequenced and the locations of the associated deep-sequencing data. Although we have technical replicates for the antibody selections for each library, we group these and just analyze variation at the level of biological replicates.

 - [./data/subamplicon_alignspecs.txt](./data/subamplicon_alignspecs.txt): the alignment specs for the [barcoded subamplicon sequencing](https://jbloomlab.github.io/dms_tools2/bcsubamp.html).
