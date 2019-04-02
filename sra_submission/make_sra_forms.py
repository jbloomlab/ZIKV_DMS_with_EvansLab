''' modified by Danny Lawrence, March 2019, to upload Zika MR766 deep mutational scanning datasets'''

import glob
import os
import sys
import pandas as pd
import tarfile

# The experiments are listed in `experiments.csv` There are six different concentrations, three biological replicates, and various technical replicates of S139/1. There are two concentrations, three biological replicates, and various technical replicates of C179. There are also four mock-selected samples (two for library 1, and one each for libraries 2 and 3), and finally a WT plasmid control.
exptsdf = pd.read_csv('./experiment_list.csv', nrows=29)
samples = exptsdf['experiment'].tolist()

# where are the FASTQ files?
origpath_dict = {'A': '/fh/fast/bloom_j/SR/ngs/illumina/bloom_lab/170609_SN367_0933_BHMLH5BCXY_lane2/Unaligned/Project_bloom_lab/',
                 'B': '/fh/fast/bloom_j/SR/ngs/illumina/bloom_lab/180427_SN367_1155_BHHJMKBCX2/Unaligned/Project_bloom_lab/',
                 'C': '/fh/fast/bloom_j/SR/ngs/illumina/bloom_lab/180928_D00300_0621_BHMJGNBCX2/Unaligned/Project_bloom_lab/'}
# the files themselves are in subdirectories of this directory, each sample name is a subdirectory.

# new biosample added (Sept 2016) to the existing bioproject of WSN HA libraries; this sample is labeled something like
# "WSN HA libraries selected with monoclonal antibodies":
biosampleid = 'SAMN05789126' # this biosampleid was created for the initial mab datasets and should be OK to use for the additional FI6v3 data.
# the bioproject number is PRJNA309339 for WSN HA sequencing.
bioprojectid = 'PRJNA309339' # this bioprojectid was created for the viruses paper when first describing the virus libraries; should be OK to continue to use this.

tar_filename = 'ZIKV_E_Ab_selection.tar'
make_tar = False # set to true to actually make the tarfile; probably want to test with *make_tar = False* at first to make sure all the files are found properly.

fnew = open('submissionform2019_ZKA64-ZKA185.tsv','w')

# header:
fnew.write('bioproject_accession\tbiosample_accession\tlibrary_ID\ttitle\tlibrary_strategy\tlibrary_source\tlibrary_selection\tlibrary_layout\tplatform\tinstrument_model\tdesign_description\tfiletype\tfilename\tfilename2\tfilename3\tfilename4\tfilename5\tfilename6\n')

# open a tar file to add all the fastq files to:
if make_tar:
    tar_out = tarfile.open(tar_filename, mode='w')

# go through samples, adding information to the submission form and adding FASTQ.gz files to the tar archive:
for sample in samples:
    print ('\nProcessing sample {0}...'.format(sample))

    FASTQ_directory = exptsdf[exptsdf['experiment'] == sample]['fastqdir'].item()
    #FASTQ_directory = origpath_dict[FASTQ_key]

    # sample_name is "library_ID" on spreadsheet
    sample_name = exptsdf[exptsdf['experiment'] == sample]['sample_name'].item()

    if ('ZKA64' in sample_name):
        title = 'MR766 E library {0} neutralized with ZKA64'.format(sample_name[-1])
    elif ('ZKA185' in sample_name):
        title = 'MR766 E library {0} neutralized with ZKA185'.format(sample_name[-1])
    elif ('Mock' in sample_name):
        title = 'MR766 E library {0} mock-neutralized'.format(sample_name[-1])
    elif ('HCV' in sample_name):
        title = 'MR766 E library {0} selected with control antibody AR3A'.format(sample_name[-1])
    elif ('zikv-' in sample_name):
        title = 'MR766 E library {0} plasmid'.format(sample_name[5:])
    elif ('zikv_' in sample_name):
        title = 'MR766 E library {0} selected in Vero cell culture'.format(sample_name[5:])


    # Each sample is a row in the spreadsheet with several filenames at the end of the row, here is everything until the filenames:
    fnew.write('{0}\t{1}\t{2}\t{3}\tAMPLICON\tVIRAL RNA\tPCR\tpaired\tILLUMINA\tIllumina HiSeq 2500\t250-nt paired-end reads of Zika virus E PCR amplicons\tfastq\t'.format(bioprojectid,biosampleid,sample_name,title))

    # Now add a tab-spaced entry in this library_ID's row for each file that will be uploaded:
    files = sorted(glob.glob('{0}Sample_{1}*/*.fastq.gz'.format(FASTQ_directory, sample_name)))
    print ('Found these files:'), files
    
    for fname in files:
        print ('Processing {0}...'.format(fname))
        sys.stdout.flush() # forces output to be written to the terminal
        short_fname = fname[fname.rfind('/')+1:]
    
        fnew.write('{0}\t'.format(short_fname))
        fnew.flush()
    
        if make_tar:
            print ('Adding {0} to the tar archive'.format(fname))
            # use the short filename when adding to the tar instead of the full path
            tar_out.add(fname, arcname=short_fname)
    
    fnew.write('\n') # end the current sample's line in the spreadsheet.

fnew.close()

if make_tar:
    tar_out.close()
    print ('Completed.')
    print ('Contents of {0}:'.format(tar_filename))
    t = tarfile.open(tar_filename, 'r')
    for member_info in t.getmembers():
        print (member_info.name)
