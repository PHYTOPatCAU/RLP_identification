import os
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess

def extract_sequences(fasta_file, protein_ids):
    """Extract sequences from a FASTA file based on protein IDs."""
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.id.split("|")[0]  # Modify if ID format differs
        if seq_id in protein_ids:
            sequences[seq_id] = record
    return sequences

def align_sequences(seq1, seq2, output_file):
    """Align two sequences using Clustal Omega."""
    with open("temp.fasta", "w") as temp_fasta:
        SeqIO.write([seq1, seq2], temp_fasta, "fasta")
    subprocess.run(["clustalo", "-i", "temp.fasta", "-o", output_file, "--force", "--dealign"], check=True)
    os.remove("temp.fasta")

def find_overlap_and_extract(seq1, seq2):
    """Find original seq in extended seq, and search upstream for stop→M→original start."""
    seq1_str = str(seq1.seq)
    seq2_str = str(seq2.seq)

    overlap_start = seq2_str.find(seq1_str)
    if overlap_start == -1:
        return seq1  # fallback if original not found in extended

    # Get upstream region
    extended_region = seq2_str[:overlap_start]

    # Search for last stop codon (*) before original start
    stop_pos = extended_region.rfind("*")
    if stop_pos == -1:
        return seq1  # No upstream stop codon

    # Search forward from stop codon for Methionine (M) before original start
    methionine_pos = extended_region.find("M", stop_pos)
    if methionine_pos != -1 and methionine_pos < overlap_start:
        extended_seq_str = seq2_str[methionine_pos:]  # from M to end
        return SeqRecord(Seq(extended_seq_str), id=seq2.id, description=seq2.description)

    return seq1  # No valid M found after stop codon and before original start

def record_gff_updates(csv_file, protein_id, extended_length, strand):
    """Record GFF coordinate updates in a CSV file."""
    with open(csv_file, "a", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([protein_id, strand, extended_length * 3])

def main():
    # Base directory
    base_dir = "~/IPS/genome.proteins/Athaliana/Correcting_the_gff_Ath/Extesnion_results"

    # Input files
    original_fasta = os.path.join(base_dir, "original_protein.fasta")
    extended_fasta = os.path.join(base_dir, "Extracting_proteins_added_with_ext_name.fasta")
    protein_ids_file = os.path.join(base_dir, "protein_ids.txt")
    output_fasta = os.path.join(base_dir, "Extracting_the_compared_and_justified_sequences.fasta")
    alignment_dir = os.path.join(base_dir, "alignments")
    csv_file = os.path.join(base_dir, "Extracting_extended_compared_justified_gff_updates.csv")

    os.makedirs(alignment_dir, exist_ok=True)

    # Prepare CSV header
    with open(csv_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Protein_ID", "Strand", "BP_To_Add"])

    # Read protein IDs
    with open(protein_ids_file, "r") as f:
        protein_ids = set(line.strip() for line in f)

    # Extract sequences
    original_sequences = extract_sequences(original_fasta, protein_ids)
    extended_sequences = extract_sequences(extended_fasta, protein_ids)

    # Process each protein ID
    with open(output_fasta, "w") as output_handle:
        for protein_id in protein_ids:
            original_seq = original_sequences.get(protein_id)
            extended_seq = extended_sequences.get(protein_id)

            if not original_seq or not extended_seq:
                continue  # Skip if sequence is missing

            # Optional: Align and save
            msa_file = os.path.join(alignment_dir, f"{protein_id}_alignment.fasta")
            align_sequences(original_seq, extended_seq, msa_file)

            # Extract extended protein
            extracted_seq = find_overlap_and_extract(original_seq, extended_seq)

            # Calculate how many extra amino acids were added
            extra_length = len(extracted_seq.seq) - len(original_seq.seq)

            # Write to output and update GFF correction file if extended
            if extra_length > 0:
                strand = "+"  # You must update this later from GFF data
                record_gff_updates(csv_file, protein_id, extra_length, strand)

            SeqIO.write(extracted_seq, output_handle, "fasta")

if __name__ == "__main__":
    main()
