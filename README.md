# Zika virus E protein deep mutational scanning
Deep mutational scanning of Zika virus E protein.
Experiments performed by the [Matt Evans lab](http://labs.icahn.mssm.edu/evanslab/).
Sequencing and computational analyses performed by Danny Lawrence and Jesse Bloom in the [Bloom lab](https://research.fhcrc.org/bloom/en.html).

## Quick summary
Look at the Jupyter notebook [analysis_notebook.ipynb](analysis_notebook.ipynb), or its markdown output at [results/summary/analysis_notebook.md](results/summary/analysis_notebook.md) for the results of the deep mutational scanning.
The growth curve data are analyzed and plotted by [plot_growth_data.ipynb](plot_growth_data.ipynb), and its markdown output is at [results/summary/plot_growth_data.md](results/summary/plot_growth_data.md).

## Running the notebooks
You can run the notebooks with the [run_nbs.bash](run_nbs.bash):

    bash run_nbs.bash

Or to submit it to the server:

    sbatch -c 16 -p largenode --mem=150000 run_nbs.bash

## Analysis and results
The analysis is performed by the Jupyter notebook [analysis_notebook.ipynb](analysis_notebook.ipynb) using [dms_tools2](https://jbloomlab.github.io/dms_tools2/).
That notebook also contains plots and descriptions of the results.

The results files are placed in the [./results](results) subdirectory.
Most of the results are not tracked in this GitHub repo, but some are.
Specifically:

  - [./results/codoncounts/](results/codoncounts) contains files that give the counts of each codon mutation for each sample from the [barcoded subamplicon sequencing](https://jbloomlab.github.io/dms_tools2/bcsubamp.html).

  - [./results/prefs/](results/prefs) has the [amino-acid preferences](https://jbloomlab.github.io/dms_tools2/prefs.html) for each library, as well as the across-library un-scaled preferences ([./results/prefs/unscaled_prefs.csv](results/prefs/unscaled_prefs.csv)) and the re-scaled preferences ([./results/prefs/rescaled_prefs.csv](results/prefs/rescaled_prefs.csv)). For most purposes, this last file is the one you want.

  - [./results/muteffects](results/muteffects) has the mutational effects calculated from the amino-acid preferences.

  - [./results/diffsel/](results/diffsel) has the [differential selection](https://jbloomlab.github.io/dms_tools2/diffsel.html) for each library, as well as the across-library average mutation differential selection ([./results/diffsel/summary_ZKA64-meanmutdiffsel.csv](results/diffsel/summary_ZKA64-meanmutdiffsel.csv) for ZKA185 and [./results/diffsel/summary_ZKA185-meanmutdiffsel.csv](results/diffsel/summary_ZKA185-meanmutdiffsel.csv) for ZKA185).

  - [./results/fracsurvive/](results/fracsurvive) has the excess (above-average) [fraction surviving](https://jbloomlab.github.io/dms_tools2/fracsurvive.html) for each library, as well as the across-library excess fraction surviving ([./results/diffsel/summary_ZKA64-meanmutfracsurvive.csv](results/diffsel/summary_ZKA64-meanmutfracsurvive.csv) for ZKA185 and [./results/diffsel/summary_ZKA185-meanmutfracsurvive.csv](results/diffsel/summary_ZKA185-meanmutfracsurvive.csv) for ZKA185).

  - [./results/logoplots](results/logoplots) contains logo plots visualizing the amino-acid preferences, differential selection, and excess fraction surviving.

  - [./results/figures](results/figures) contains plots that we anticipate being figures in the paper.

## Input data
The input data are in the [./data/](data) subdirectory.
These data consist of:

 - [./data/E.fasta](data/E.fasta): coding sequence of E protein from ZIKV MR766 strain.

 - [./data/samplelist.csv](data/samplelist.csv): all the samples that we sequenced and the locations of the associated deep-sequencing data. Although we have technical replicates for the antibody selections for each library, we group these and just analyze variation at the level of biological replicates.

 - [./data/subamplicon_alignspecs.txt](./data/subamplicon_alignspecs.txt): the alignment specs for the [barcoded subamplicon sequencing](https://jbloomlab.github.io/dms_tools2/bcsubamp.html).

 - [./data/E_alignment.fasta](data/E_alignment.fasta): alignment of ZIKV E protein created by Danny Lawrence. Downloaded from [NCBI Virus Variation Resource](https://www.ncbi.nlm.nih.gov/genome/viruses/variation/Zika/) all ZIKV sequences (2019-03-26) with the query settings as follows:
   - Sequence type: Protein
   - Host: any
   - Region/Country: any
   - Genome region: E
   - Isolation source: any
   - Collapse identical sequences: checked ("Note: All groups of identical sequences in the dataset will be represented by the oldest sequence in the group.")
   - Added the E region of our WT sequence for [MR766](https://www.ncbi.nlm.nih.gov/nuccore/KX830961) as the reference sequence and used [`phydms_prepalignment`](http://jbloomlab.github.io/phydms/phydms_prepalignment.html) to trim the alignment and remove any sequences that are redundant or incomplete. A number of sequences appeared to be derived from an experiment or otherwise spurious and were removed manually. 

 - [data/domains.csv](data/domains.csv): the domain structure of E as defined in Figure 1 of [Dai et al (2016)](https://www.cell.com/cell-host-microbe/fulltext/S1931-3128(16)30149-4).

 - [./data/all_growth_data.csv](data/all_growth_data.csv):
   Growth curve data generated by Matt Evans.

 - PDB structures of the E protein:

   - [data/5ire.pdb](data/5ire.pdb):
     PDB [5ire](http://www.rcsb.org/structure/5ire)

   - [data/5ire_monomer.pdb](data/5ire_monomer.pdb):
     PDB [5ire](http://www.rcsb.org/structure/5ire), chain A only

   - [data/5u4w.pdb](data/5u4w.pdb):
     PDB [5u4w](http://www.rcsb.org/structure/5u4w)

   - [data/6co8.pdb](data/6co8.pdb):
     PDB [5co8](http://www.rcsb.org/structure/6co8)

   - [data/6co8_monomer.pdb](data/6co8_monomer.pdb):
     PDB [5co8](http://www.rcsb.org/structure/6co8), chain A only

 - Secondary structure and solvent accessibility calculations made by running the [dssp webserver](http://www.cmbi.ru.nl/xssp/) on the above structures:

   - [data/5ire.dssp](data/5ire.dssp)
   
   - [data/5ire_monomer.dssp](data/5ire_monomer.dssp)
   
   - [data/6co8.dssp](data/6co8.dssp)
   
   - [data/6co8_monomer.dssp](data/6co8_monomer.dssp)

 - Plasmid maps:

   - [data/1725_ZIKV_MR766_int_WT.gb](data/1725_ZIKV_MR766_int_WT.gb): plasmid for wild-type ZIKV reverse genetics.

   - [data/1726_ZIKV_MR766_int_GFP.gb](data/1726_ZIKV_MR766_int_GFP.gb): recipient plasmid for cloning mutant libraries.

## PyMol scripts
The subdirectory [pymol_scripts](pymol_scripts) contains scripts for generating PyMol structure images.
