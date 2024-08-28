import os
import random
from pymol import cmd

# List of standard amino acids
amino_acids = [
    'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
    'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'
]

# Function to generate a random peptide sequence
def generate_peptide_sequence(length=13):
    return ''.join(random.choices(amino_acids, k=length))

# Function to generate a peptide in PyMOL and save as PDB
def generate_peptide_pdb(sequence, filename):
    # Create a peptide using PyMOL's built-in functions
    cmd.fab(sequence, "peptide", ss=1)
    
    # Save the peptide as a PDB file
    cmd.save(filename, "peptide")
    cmd.delete("all")

# Main function to generate 100 peptides and save them in different folders
def generate_and_save_peptides(num_peptides=100, sequence_length=13):
    for i in range(num_peptides):
        sequence = generate_peptide_sequence(sequence_length)
        directory = f"peptide_{i + 1}"  # Each sequence in its own folder
        os.makedirs(directory, exist_ok=True)
        pdb_filename = os.path.join(directory, f"{sequence}.pdb")
        fasta_filename = os.path.join(directory, f"{sequence}.fasta")
        
        # Generate the PDB file using PyMOL
        generate_peptide_pdb(sequence, pdb_filename)
        
        # Save the sequence as a FASTA file
        with open(fasta_filename, "w") as fasta_file:
            fasta_file.write(f">{sequence}\n")
            fasta_file.write(sequence)
        
        print(f"Saved PDB and FASTA for {sequence} in {directory}")

# Generate 100 peptides and save them
generate_and_save_peptides(100, 13)

