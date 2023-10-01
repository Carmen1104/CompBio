# Define a dictionary mapping DNA codons to amino acids
codon_table = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}

# Function to translate DNA to mRNA and then to amino acids
def translate_dna_to_amino_acids(dna_sequence):
    # Replace "T" with "U" to get mRNA sequence
    mrna_sequence = dna_sequence.replace("T", "U")
    
    # Initialize an empty string to store the amino acid sequence
    amino_acid_sequence = ""
    
    # Initialize a dictionary to store codon frequencies
    codon_frequencies = {}
    
    # Iterate through the mRNA sequence in triplets (codons)
    for i in range(0, len(mrna_sequence) - 2, 3):
        codon = mrna_sequence[i:i + 3]
        amino_acid = codon_table.get(codon, "X")  # "X" for unknown codons
        amino_acid_sequence += amino_acid
        
        # Count the frequency of each RNA codon
        codon_frequencies[codon] = codon_frequencies.get(codon, 0) + 1
    
    return amino_acid_sequence, codon_frequencies

# Input DNA sequence
input_dna = "AAUGCUAAU"

# Input Aminoacid
input_aminoacid = "NAN"

# Translate DNA to mRNA and then to amino acids
amino_acid_sequence, codon_frequencies = translate_dna_to_amino_acids(input_dna)

# Output the results
print("Input Aminoacid =", input_aminoacid)
for codon, frequency in codon_frequencies.items():
    print(f"{codon} = {frequency}")
