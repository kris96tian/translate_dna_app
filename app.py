from flask import Flask, render_template, request

app = Flask(__name__)


def translate_dna_to_protein(dna_sequence):
    genetic_code = {
        'TTT': ('F', 'Phenylalanine'), 'TTC': ('F', 'Phenylalanine'),
        'TTA': ('L', 'Leucine'), 'TTG': ('L', 'Leucine'),
        'TCT': ('S', 'Serine'), 'TCC': ('S', 'Serine'), 'TCA': ('S', 'Serine'), 'TCG': ('S', 'Serine'),
        'TAT': ('Y', 'Tyrosine'), 'TAC': ('Y', 'Tyrosine'),
        'TGT': ('C', 'Cysteine'), 'TGC': ('C', 'Cysteine'), 'TGG': ('W', 'Tryptophan'),
        'CTT': ('L', 'Leucine'), 'CTC': ('L', 'Leucine'), 'CTA': ('L', 'Leucine'), 'CTG': ('L', 'Leucine'),
        'CCT': ('P', 'Proline'), 'CCC': ('P', 'Proline'), 'CCA': ('P', 'Proline'), 'CCG': ('P', 'Proline'),
        'CAT': ('H', 'Histidine'), 'CAC': ('H', 'Histidine'),
        'CAA': ('Q', 'Glutamine'), 'CAG': ('Q', 'Glutamine'),
        'CGT': ('R', 'Arginine'), 'CGC': ('R', 'Arginine'), 'CGA': ('R', 'Arginine'), 'CGG': ('R', 'Arginine'),
        'ATT': ('I', 'Isoleucine'), 'ATC': ('I', 'Isoleucine'), 'ATA': ('I', 'Isoleucine'),
        'ATG': ('M', 'Methionine'),
        'ACT': ('T', 'Threonine'), 'ACC': ('T', 'Threonine'), 'ACA': ('T', 'Threonine'), 'ACG': ('T', 'Threonine'),
        'AAT': ('N', 'Asparagine'), 'AAC': ('N', 'Asparagine'),
        'AAA': ('K', 'Lysine'), 'AAG': ('K', 'Lysine'),
        'AGT': ('S', 'Serine'), 'AGC': ('S', 'Serine'),
        'AGA': ('R', 'Arginine'), 'AGG': ('R', 'Arginine'),
        'GTT': ('V', 'Valine'), 'GTC': ('V', 'Valine'), 'GTA': ('V', 'Valine'), 'GTG': ('V', 'Valine'),
        'GCT': ('A', 'Alanine'), 'GCC': ('A', 'Alanine'), 'GCA': ('A', 'Alanine'), 'GCG': ('A', 'Alanine'),
        'GAT': ('D', 'Aspartic Acid'), 'GAC': ('D', 'Aspartic Acid'),
        'GAA': ('E', 'Glutamic Acid'), 'GAG': ('E', 'Glutamic Acid'),
        'GGT': ('G', 'Glycine'), 'GGC': ('G', 'Glycine'), 'GGA': ('G', 'Glycine'), 'GGG': ('G', 'Glycine'),
        'TAA': ('*', 'Stop'), 'TAG': ('*', 'Stop'), 'TGA': ('*', 'Stop')
    }

    dna_sequence = dna_sequence.upper().replace(' ', '')  # Ensure DNA sequence is uppercase and remove spaces
    if not all(nucleotide in 'ATCG' for nucleotide in dna_sequence):
        raise ValueError("Invalid DNA sequence. Only 'A', 'T', 'C', and 'G' are allowed.")

    protein_sequence = []
    amino_acid_sequence = []

    for i in range(0, len(dna_sequence), 3):
        codon = dna_sequence[i:i + 3]
        if len(codon) < 3:
            break  # Skip incomplete codons at the end of the sequence

        if codon in genetic_code:
            amino_acid, aa_name = genetic_code[codon]
            protein_sequence.append(amino_acid)
            amino_acid_sequence.append(aa_name)
        else:
            amino_acid_sequence.append('Unknown')

    protein_sequence = ''.join(protein_sequence)
    amino_acid_sequence = '-'.join(amino_acid_sequence)

    return protein_sequence, amino_acid_sequence


@app.route('/', methods=['GET', 'POST'])
def index():
    protein_sequence = ''
    amino_acid_sequence = ''
    error_message = ''

    if request.method == 'POST':
        dna_sequence = request.form['dna_sequence']
        try:
            protein_sequence, amino_acid_sequence = translate_dna_to_protein(dna_sequence)
        except ValueError as e:
            error_message = str(e)

    return render_template('index.html', protein_sequence=protein_sequence, amino_acid_sequence=amino_acid_sequence, error_message=error_message)

if __name__ == '__main__':
    app.run(debug=True)