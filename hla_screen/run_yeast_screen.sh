#!/bin/bash
while getopts "i:o:d:" flag; do
    case "${flag}" in
        i) input="${OPTARG}" ;;
        o) output="${OPTARG}" ;;
        d) directory="${OPTARG}" ;;
    esac
done
echo "input: $input"
echo "output: $output"
echo "directory: $directory"

# Make the working directory the one for the current script:
cd "$(dirname "$0")"
echo "Current working directory: $(pwd)"

# Read config file:
. ./config.sh

# Process Fastq file:
sh ./process_fastq_pipeline.sh -i $input -o $output -d $directory

# Decode ORFs and match with the reference sequence
python ./translate_sequences_pipeline.py $directory/$input'_extracted_consensus.fasta' $directory/$input'_extracted_consensus_proteins.fasta' ./reference/$ref_prot_fasta