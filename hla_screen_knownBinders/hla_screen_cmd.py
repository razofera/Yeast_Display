import pandas as pd
import configparser

def extract_sequences_from_fastq(file_path, start_char, stop_char):
    sequences = []
    with open(file_path, 'r') as file:
        while True:
            header = file.readline().strip()
            if not header:
                break
            sequence = file.readline().strip()
            plus_line = file.readline().strip()
            quality = file.readline().strip()
            
            # Extract the subsequence from start_char to stop_char
            subsequence = sequence[start_char:stop_char]
            sequences.append(subsequence)
    
    return sequences

def count_sequence_occurrences(sequence_list):
    sequence_counts = {}
    for sequence in sequence_list:
        if sequence in sequence_counts:
            sequence_counts[sequence] += 1
        else:
            sequence_counts[sequence] = 1
    return sequence_counts

def extract_keys(dictionary):
    keys = list(dictionary.keys())
    return keys

def read_fasta(file_path):
    sequences = {}
    with open(file_path, 'r') as file:
        header = None
        sequence = ''
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if header:
                    sequences[header] = sequence
                header = line[1:]  # Remove '>'
                sequence = ''
            else:
                sequence += line
        if header:
            sequences[header] = sequence
    return sequences

def find_matching_sequences(fasta_sequences, query_sequences):
    matched_sequences = {}
    for header, sequence in fasta_sequences.items():
        for query in query_sequences:
            if query in sequence:
                matched_sequences[sequence] = header
                break
    return matched_sequences

def extract_unique_values(dictionary):
    unique_values = set(dictionary.values())
    return list(unique_values)
    
def find_matching_header(fasta_sequences, query_sequences):
    matched_sequences = {}
    for header, sequence in fasta_sequences.items():
        for query in query_sequences:
            if query in header:
                matched_sequences[query] = sequence
                break
    return matched_sequences

def extract_subsequences(dictionary, start_char, stop_char):
    modified_dict = {}
    for key, value in dictionary.items():
        subsequence = key[start_char:stop_char]
        modified_dict[subsequence] = value
    return modified_dict

config = configparser.ConfigParser()
config.read('./config.ini')

# Filepaths:
file_path = config['DEFAULT']['file_path']
fasta_file_path = config['DEFAULT']['fasta_file_path']
peptide_file_path = config['DEFAULT']['peptide_file_path']
csv_file_path = config['DEFAULT']['csv_file_path']

# Variables: Adjust these depending on how long adaptors and barcodes are in sequences and in reference file
start_char = config['DEFAULT']['start_char']
stop_char = config['DEFAULT']['stop_char']
start_ref = config['DEFAULT']['start_ref']
stop_ref = config['DEFAULT']['stop_ref']

def main():
    # Extract sequences and remove adaptors in fastq files:
    sequences = extract_sequences_from_fastq(file_path, start_char, stop_char)
    
    # Count number of occurances of each sequence from screen:
    sequence_counts = count_sequence_occurrences(sequences)
    
    # Check if sequences match reference list:
    fasta_sequences = read_fasta(fasta_file_path)
    query_sequences = extract_keys(sequence_counts)
    matched_sequences = extract_subsequences(
        find_matching_sequences(fasta_sequences, query_sequences),
        start_ref,
        stop_ref)
    
    # Check the peptide sequence for each query sequence:
    query_headers = extract_unique_values(matched_sequences)
    peptide_sequences = read_fasta(peptide_file_path)
    matched_peptides = find_matching_header(peptide_sequences, query_headers)
    
    # Convert dictionaries to DataFrames and merge:
    df1 = pd.DataFrame(list(sequence_counts.items()), columns=['Sequence', 'Count'])
    df2 = pd.DataFrame(list(matched_sequences.items()), columns=['Sequence', 'Identifier'])
    df3 = pd.DataFrame(list(matched_peptides.items()), columns=['Identifier', 'Protein'])
    merged_df = pd.merge(df1, df2, on='Sequence')
    final_df = pd.merge(merged_df, df3, on='Identifier')
    
    # Save DataFrame to CSV file
    final_df.to_csv(csv_file_path)
    print(f"Results have been saved to {csv_file_path}")
    
if __name__ == "__main__":
    main()