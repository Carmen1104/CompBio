# Define a dictionary mapping DNA codons to amino acids
codon_table = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}

# Function to get the complement of a DNA sequence
def get_complement(dna_sequence):
    complement_dict = {"A": "T", "T": "A", "C": "G", "G": "C"}
    complement_sequence = "".join(complement_dict[base] for base in dna_sequence)
    return complement_sequence

# Function to translate DNA to mRNA and then to amino acids
def translate_dna_to_protein(dna_sequence):
    # Replace "T" with "U" to get mRNA sequence
    mrna_sequence = dna_sequence.replace("T", "U")
    
    # Initialize an empty string to store the amino acid sequence
    protein_sequence = ""
    
    # Iterate through the mRNA sequence in triplets (codons)
    for i in range(0, len(mrna_sequence), 3):
        codon = mrna_sequence[i:i + 3]
        amino_acid = codon_table.get(codon, "X")  # "X" for unknown codons
        protein_sequence += amino_acid
    
    return protein_sequence

# Input DNA sequence
input_dna = "TTACGA"

# Calculate the complement sequence
complement_dna = get_complement(input_dna)

# Translate DNA to mRNA and then to amino acids
mRNA_sequence = input_dna.replace("T", "U")
amino_acid_sequence = translate_dna_to_protein(input_dna)

# Output the results
print("Input DNA =", input_dna)
print("Complement =", complement_dna)
print("mRNA =", mRNA_sequence)
print("Aminoacid =", amino_acid_sequence)
