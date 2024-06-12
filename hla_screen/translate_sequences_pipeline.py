from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW
import os
import glob
from Bio.Seq import Seq
import sys
input_file = sys.argv[1]
output_file = sys.argv[2]
ref_fasta = sys.argv[3]

print(input_file)
print(output_file)
print(ref_fasta)

def translate_seq_to_ORFs(input_record):
    s = input_record.seq
    # Generate all translations from 6 ORFs possible (3 forward, 3 reverse)
    out = []
    for i in range(3):
        aa_seq = s[i:].translate(to_stop=True)
        aa_seq_rev = s.reverse_complement()[i:].translate(to_stop=True)
        out.append(aa_seq)
        out.append(aa_seq_rev)
    return(out)

def extract_id_seq(accession_id, reference_fasta):
    # Read the accession IDs from reference fasta
    with open(reference_fasta, 'r') as f2:
        records = list(SeqIO.parse(f2, 'fasta'))
        accession_id = accession_id.split('|')[0]
        for record in records:
            if record.id.split('|')[0] == accession_id:
                target_sequence = record.seq
                break
        f2.close()
        return(target_sequence)

def match_seqs(seq_list, ref_seq):
    matched_string = []
    # Iterate over each string in the list
    for query in seq_list:
        # Check if the substring is present in the string
        if query in ref_seq:
            # Append the matching string to the list
            matched_string.append(query)
    if len(matched_string):
        matched_string = matched_string[0]
    return matched_string

# Usage:

for record in SeqIO.parse(input_file, "fasta"):
    orfs = translate_seq_to_ORFs(record)
    seq_id = record.id
    ref_seq = extract_id_seq(seq_id, ref_fasta)
    match_seq = match_seqs(orfs, ref_seq)
    # Store the protein sequence in a new FASTA record
    protein_record = SeqRecord(match_seq, id=record.id)
    SeqIO.write(protein_record, output_file, "fasta")
