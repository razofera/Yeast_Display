# Yeast_Display
Tools for Yeast Display Screen Binding Analysis

READ ME

# INSTALLATION
- clone repository from github
- optional: add 'run_yeast_screen.sh' to the PATH variable
- cd to the download directory. Run 'sh setup.sh' to setup conda environment
- change the variables in the 'config.sh' file to your specific system
	this includes putting your own reference file fasta into the 'reference' folder and updating the variable with the name

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

""" Example Code """
cd ./path/to/run_yeast_screen.sh/folder /
conda activate yeast_screen /
sh run_yeast_screen.sh -i "test_input_2peptides" -o "output" -d "/mnt/c/Users/razof/Documents/linux/aly/data/hla_screen/tests"

"""              """
