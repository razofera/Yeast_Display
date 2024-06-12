#!/bin/bash
while getopts "i:o:d:" flag; do
    case "${flag}" in
        i) input="${OPTARG}" ;;
        o) output="${OPTARG}" ;;
        d) directory="${OPTARG}" ;;
    esac
done
# echo "input: $input"
# echo "output: $output"
# echo "directory: $directory"

# Read config file:
. ./config.sh

# Step 0: Run FastQC on all FASTQ files (you can use *.fq wildcard)
fastqc $directory/$input'.fastq' -o $directory

# Step 0.1: Run Trimgalore (https://docs.tinybio.cloud/docs/trim-galore-tutorial)
trim_galore $directory/$input'.fastq' -o $directory

# Step 1: Align reads
# Check first if the reference has already been indexed
if [ -f "./reference/"$ref_seq_fasta".pac" ]; then
    echo "File exists. Skipping the code."
else
    echo "File does not exist. Proceeding with the code."
    ./bwa-mem2-2.2.1_x64-linux/bwa-mem2 index ./reference/$ref_seq_fasta
fi
./bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem ./reference/$ref_seq_fasta $directory/$input'_trimmed.fq' > $directory/$input'.sam'

# Step 2: Convert SAM to sorted BAM
samtools view -Sb $directory/$input'.sam' | samtools sort -o $directory/$input'_sorted.bam'

# Step 3: Call peaks with MACS2
macs2 callpeak -t $directory/$input'_sorted.bam' -f BAM -g hs -n peaks --nomodel --extsize 147 --outdir $directory

# Step 4: Generate consensus sequence
samtools view -b -h -L $directory/peaks_peaks.narrowPeak $directory/$input'_sorted.bam' > $directory/$input'_sam_regions.bam'
samtools mpileup -uf ./reference/$ref_seq_fasta $directory/$input'_sam_regions.bam' | bcftools call -m -O z - > $directory/$input'.vcf.gz'
bcftools index $directory/$input'.vcf.gz'
bcftools consensus -f ./reference/$ref_seq_fasta $directory/$input'.vcf.gz' > $directory/$input'_consensus.fasta'

# Extract the consensus sequences
bedtools getfasta -fi $directory/$input'_consensus.fasta' -bed $directory/peaks_peaks.narrowPeak -fo $directory/$input'_extracted_consensus.fasta'