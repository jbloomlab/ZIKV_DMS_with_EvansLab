''' modified by Danny Lawrence, March 2019, to upload Zika MR766 deep mutational scanning datasets'''

import glob
import os
import sys
import pandas as pd
import tarfile

# The experiments are listed in `samplelist.csv` 
exptsdf = (pd.read_csv('../data/samplelist.csv')
			.assign(sample_name=lambda x: x.selection + '_' + x.library))


biosampleid = 'SAMN11335355' # this biosampleid was created for the initial ZIKV DMS and mab selection.
# the bioproject number is PRJNA530795 for ZIKV DMS sequencing.
bioprojectid = 'PRJNA530795' 

tar_filename = 'ZIKV_E_Ab_selection.tar'
make_tar = True # set to true to actually make the tarfile; probably want to test with *make_tar = False* at first to make sure all the files are found properly.

fnew = open('submissionform2019_ZKA64-ZKA185.tsv','w')

# header:
fnew.write('bioproject_accession\tbiosample_accession\tlibrary_ID\ttitle\tlibrary_strategy\tlibrary_source\tlibrary_selection\tlibrary_layout\tplatform\tinstrument_model\tdesign_description\tfiletype\tfilename\tfilename2\tfilename3\tfilename4\tfilename5\tfilename6\tfilename7\tfilename8\tfilename9\tfilename10\tfilename11\tfilename12\tfilename13\tfilename14\tfilename15\tfilename16\tfilename17\tfilename18\tfilename19\tfilename20\tfilename21\tfilename22\tfilename23\tfilename24\n')

# open a tar file to add all the fastq files to:
if make_tar:
    tar_out = tarfile.open(tar_filename, mode='w')

# go through samples, adding information to the submission form and adding FASTQ.gz files to the tar archive:
for sample_tup in exptsdf.itertuples():
    sample_name = sample_tup.sample_name
    print ('\nProcessing sample {0}...'.format(sample_name))

    FASTQ_directory = sample_tup.R1
    print(FASTQ_directory)

    # sample_name is "library_ID" on spreadsheet
    selection = sample_tup.selection
    library = sample_tup.library

    if ('ZKA64' in sample_name):
        title = 'MR766 E library {0} neutralized with ZKA64'.format(library)
    elif ('ZKA185' in sample_name):
        title = 'MR766 E library {0} neutralized with ZKA185'.format(library)
    elif ('no-antibody' in sample_name):
        title = 'MR766 E library {0} mock-neutralized'.format(library)
    elif ('control-antibody' in sample_name):
        title = 'MR766 E library {0} selected with control antibody AR3A'.format(library)
    elif ('plasmid' in sample_name):
        title = 'MR766 E library {0} plasmid'.format(sample_name[5:])
    elif ('virus' in sample_name):
        title = 'MR766 E library {0} selected in Vero cell culture'.format(library)


    # Each sample is a row in the spreadsheet with several filenames at the end of the row, here is everything until the filenames:
    fnew.write('{0}\t{1}\t{2}\t{3}\tAMPLICON\tVIRAL RNA\tPCR\tpaired\tILLUMINA\tIllumina HiSeq 2500\t250-nt paired-end reads of Zika virus E PCR amplicons\tfastq\t'.format(bioprojectid,biosampleid,sample_name,title))

    # Now add a tab-spaced entry in this library_ID's row for each file that will be uploaded:
    files = sorted(glob.glob(FASTQ_directory)+glob.glob(FASTQ_directory.replace('R1', 'R2')))
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
