# Notes for uploading deep sequencing files to the NCBI SRAupload

1. Log in to the [SRA Submission Portal Wizard](https://submit-ncbi-nlm-nih-gov.offcampus.lib.washington.edu/subs/sra/) and start a new submission to acquire new BioProject and BioSample numbers for the project. Once these numbers are given, edit the script `make_sra_forms.py` to write this information into the metadata submission form.
2. Run `make_sra_forms.py` to write the metadata submission form (.tsv) and to create the .tar file to upload to the SRA.
3. Pre-upload sequencing files using FTP.
    - Navigate to the directory containing the .tar file
    - Type the following command to establish an FTP connection: `ftp ftp-private.ncbi.nlm.nih.gov`;
      username: `subftp`
      password:  `w4pYB9VQ`
    - Navigate to your account folder (e.g., `cd uploads/djpl_uw.edu_1FRr19Sn`)
    - Create a new subfolder by typing `mkdir ZIKV_sequencing` and navigate to this new folder (`cd ZIKV_sequencing`)
    - Copy tar file by typing `put ZIKV_E_Ab_selection.tar`
4. Once the preupload is complete, return to the [SRA Submission Portal Wizard](https://submit.ncbi.nlm.nih.gov/subs/sra/) and start a new submission using the biosample and bioproject numbers assigned previously.
5. Write in personal and general information for the submission and upload the metadata .tsv file
6. In the "Files" tab of the submission, click "Select preload folder" to select the .tar file that has just been uploaded
