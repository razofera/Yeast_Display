# Yeast_Display
Tools for Yeast Display Screen Binding Analysis

READ ME

# KNOWN REFERENCE
- If a reference file with sequences being queried in the Yeast Display Screen is available, you can run the 'hla_screen_knownBinders' pipeline.

  # INSTALLATION
  - clone 'hla_screen_knownBinders' directory
  
  # DEPENDENCIES
  - Python 3
  - jupyter notebook
  - Pandas Library
 
  # TESTING
  - Test the code by running it on the provided test files
  - 	test_known_unlabeled - fastq file with sequencing result
  - 	test_known_binders - fasta file with reference seqeunces from yeast display screen
  - 	test_known_peptides - fasta file with reference peptide sequences
  - The output from the code should match:
  - 	test_merged_df.csv
 
  # RUNNING
  - Adjust the variables for filepaths as well as the 'start' and 'stop' variables for adaptors/barcodes in your sequencing library

_________________________________________________________________________________________________________________________________________________________________

# UNKNOWN REFERENCE
- If there are no reference files, then you can run the following pipeline
- The pipeline first constructs consensus sequences for the peptide binders in the screen, then translates possible ORFs and finally matches to known peptides

# INSTALLATION
- clone repository from github
- cd to the download directory 'hla_screen' folder. Run 'sh setup.sh' to setup conda environment
- download reference sequence fasta files (CCDS_nucleotide.current.fna.gz, CCDS_protein.current.faa.gz) from https://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/
- 	unzip them using 'gunzip CCDS........fasta.gz'
- 	rename them to CCDS.fasta and CCDS_protein.fasta, respectively and put them inside a new folder called 'reference' inside the 'hla_screen' directory
- change the variables in the 'config.sh' file to your specific system, including the names of the CCDS fasta files

# DEPENDENCIES
- OS: Ubuntu 22.04 LTS  
	fastqc (v0.11.9)  
	trim_galore (v0.6.10)  
	samtools (v1.13) (using htslib 1.13+ds)  
	macs2 (v2.2.9.1)  
	bcftools (v1.13)  
	bedtools (v2.30.0)  
- Python: Anaconda3 (v24.1.2)  
	Biopython (v1.83)  

# RUNNING THE PIPELINE
The pipeline is run using the script 'run_yeast_screen.sh'
- cd to the folder that contains the script
- activate the 'yeast_screen' conda environment
- run the script with the following options:
	-i	input data filename
	-o	output data filename
	-d	directory where input/output files are

Example Test Code:  
cd /path/to/run_yeast_screen.sh/folder  
sh setup.sh  
conda activate yeast_screen  
sh run_yeast_screen.sh -i "insulin_49mer_25bpoverlap.fastq" -o "output" -d "./test"  
  
The output of this code should match the 'test_output.fasta' file

