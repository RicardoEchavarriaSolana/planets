
import biotite 
import biotite.sequence as seq 
import biotite.sequence.graphics as graphics
import matplotlib 
import matplotlib.pyplot as plt 
matplotlib.use('Qt5Agg')
plt.ion()
plt.rcParams['font.size'] = 20
plt.rcParams['lines.linewidth'] = 2

#%% Sequence manipulation 

dir(biotite)
dir(seq)

# General sequence
dna = seq.NucleotideSequence('ATCTGCAACTGA')


# Get all the letters present in a sequence
dna.alphabet

# Get complement of a sequence 

print('\n\nThe complement of', dna, 'is', dna.complement())

# Create a protein sequence 
protein_seq = seq.ProteinSequence("MIT")
print("\n\nTranslation of", dna, "gives:", protein_seq)
protein_seq.alphabet # lists all possible symbols
protein_seq.get_symbol_frequency() # frequencies

# Convert aa's from 1-letter to 3-letters
print('\nConvert aa symbols from 1 to 3 letters: ')
for aa in protein_seq:
    print(aa, seq.ProteinSequence.convert_letter_1to3(aa))


# 1. Find the first position of Base or Nucleotide
# NOTE: .find, .rfind, .count functions
# require for the sequence to be
# in string format 

# Search from the left 
str(dna).find('A')

# Search form the right
str(dna).rfind('A')

# Count freq of nt
str(dna).count('A')

# 2. Get total nt count
dna.get_symbol_frequency() 

# 3. Search for subseq
sub_seq = seq.NucleotideSequence('TGC')
seq.find_subsequence(dna, sub_seq)

# 4. Find the ocurrence of nt in seq
print("\n\nAdenine nucleotides occur at indexes:", seq.find_symbol(dna,"A"), "of sequence", dna)
#print(seq.find_symbol_last(dna,"A"))

# 5. Plot nt count
dna_freq = dna.get_symbol_frequency() # type dictionary
plt.bar(dna_freq.keys(), dna_freq.values())
plt.pause(0.5)
plt.show(block=True)

#%% Transcription and Translation

protein = dna.translate(complete=True)

# aa's CodonTable

#print(seq.CodonTable.load(1))
table= seq.CodonTable.default_table() 

codon1= 'GAA'
aa1  = table['GAA']

print('\nCodon',codon1, 'codes for:', aa1)

# You can do the reverse too 
table['E']

# Check for Codon Table for other species
bact_ctab = seq.CodonTable.load(11)
yeast_mt_ctab = seq.CodonTable.load("Yeast Mitochondrial")

#%% Sequence allignment 

# 1. Global: Needleman-Wunsch algorithm

# 2. Local: Smith-Waterman algorithm; often preferable

import biotite.sequence.align as align

dir(align)

# Create Standard subMatrix
n_matrix = align.SubstitutionMatrix.std_nucleotide_matrix()
p_matrix = align.SubstitutionMatrix.std_protein_matrix() 

# Load subMatrix for internal database 
alpha = seq.ProteinSequence.alphabet
matrix_bl50 = align.SubstitutionMatrix(alpha,alpha, "BLOSUM50")

#Prot seq alignment

pseq1 = seq.ProteinSequence("MITTITE")
pseq2 = seq.ProteinSequence("ITTITE")

# Local alignments
loc_alignments = align.align_optimal(pseq1,pseq2,p_matrix,local=True)
print('\n\nLocal alignment for protein sequences {seq1} and {seq2}: \n'.format(seq1=pseq1, seq2=pseq2))
for a in loc_alignments:
    print(a)

fig, ax = plt.subplots(figsize=(2.0, 0.8))
graphics.plot_alignment_similarity_based(
    ax, loc_alignments[0], matrix=p_matrix, symbols_per_line=len(loc_alignments[0])
) 
fig.tight_layout() 
plt.pause(0.5)
plt.show()

#score 
alignments = loc_alignments[0]
print('\nScore:', alignments.score)
print('Recalc score:', align.score(alignments,matrix=p_matrix))
print('Seq. identity:', align.get_sequence_identity(alignments))


# Global alignments
glob_alignments = align.align_optimal(pseq1,pseq2, p_matrix,local=False)
print('\n\nGlobal alignment for protein sequences {seq1} and {seq2}: \n'.format(seq1=pseq1, seq2=pseq2))
for a in glob_alignments:
    print(a)

fig, ax = plt.subplots(figsize=(2.0, 0.8))
graphics.plot_alignment_similarity_based(
    ax,glob_alignments[0], matrix=p_matrix, symbols_per_line=len(glob_alignments[0])
) 
fig.tight_layout() 
plt.pause(0.5)
plt.show()

#score 
alignments = glob_alignments[0]
print('\nScore:', alignments.score)
print('Recalc score:', align.score(alignments,matrix=p_matrix))
print('Seq. identity:', align.get_sequence_identity(alignments))

import numpy as np 

# This is a test for github 

