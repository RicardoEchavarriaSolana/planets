
# Load pkgs

import Bio 
from IPython import get_ipython
get_ipython().magic('reset -sf')

# Note to self: 'Alphabet' feature is obsolete in current version
dir(Bio)

# Classes: 'Data' 'Seq'

from Bio.Seq import Seq 

# Seq Class methods (Functions stored in) and attributes (Variables stored in)
dir(Seq)

# Create a simple sequence

seq_1= Seq('ATGATCTCGTAA')
protein_a = Seq('MIT')

# Sequence manipulation 

# 1. Indexing/Slicing:  
# length of sequence
len(seq_1)

# Slicing of a sequence
seq_1[0:3]

# Reverse
seq_1[::-1]

# 2. Join 2 Sequences:
seq_2 = Seq('ATGATCTCGTGG')
seq_3 = seq_1[0:6]+seq_2
seq_3 += 'TAG' # add a stop codon

# 3. Find position of nucleotide in a Sequence
# first from the left
seq_3.find('ATG', start=5, end=9)
seq_3.index('ATG')

# first from the right
seq_3.rfind('G')
seq_3.rindex('G')

# All existing positions (custom-made, see nt_search() instead)
def all_pos(seq, string, n=0):
    
    if len(string) == 1:
        pos = [idx+n for idx, item in enumerate(seq) if string in item]
    else:
        pos = []
        tmp_indx = -len(string)
        for i in range(seq.count(string)):
            tmp_indx = seq.find(string, start=tmp_indx+len(string), end=len(seq))
            pos.append(tmp_indx+n)
    return pos 

all_pos(seq_3, 'TC')
# This function will return all positions for a given nt, or all the positions
#  of the first nt of a given string, with n-based indexing (defaults to 0)


# 4. Count nucleotide number
seq_3.count('G')
seq_3.count('ATG')

# 5. Find the number of repeats of a given sequence 
### Plot Frequencies ###
import matplotlib 
import matplotlib.pyplot as plt 
from collections import Counter

matplotlib.use('Qt5Agg')
plt.ion()
plt.rcParams['font.size'] = 20
plt.rcParams['lines.linewidth'] = 2

# All nts frequencies 
nt_freq = Counter(seq_3)

### Repeats of a given codon ### 
codon_list = ['TC', 'ATG']
codon_freq = []
for n in range(len(codon_list)):
    tmp_freq = seq_3.count(codon_list[n]) 
    codon_freq.append(tmp_freq)
    
# plot 
plt.bar(nt_freq.keys(),  nt_freq.values())
plt.bar(codon_list, codon_freq)
plt.ylabel('count')
plt.xlabel('seq={seq}'.format(seq=seq_3))

plt.pause(0.5)
plt.show(block=True)

#%% Protein synthesis

# template strand (complemetary to coding strand)
tseq_3 = seq_3.complement()

# Transcription (coding strand as input) 
mrna_seq_3 = seq_3.transcribe() 

# Back transcription
back_seq = mrna_seq_3.back_transcribe()
back_seq == seq_3

# Translation: 1 (single letter aas)
mrna_seq_3.translate() 

# Translation: 2 (single letter aas)
aa = seq_3.translate() 

# Create custom stop codon symbol
mrna_seq_3.translate(stop_symbol='x')

# Translation: 3 (3-letter aa's system)
from Bio.SeqUtils import seq3, seq1
seq3(aa)
seq1(seq3(aa))

# View our codonTable
from Bio.Data import CodonTable

#methods 
dir(CodonTable)

print(CodonTable.unambiguous_dna_by_name['Standard'])
print(CodonTable.unambiguous_rna_by_name['Standard'])

#%% DNA Composition: GC and AT Content

from Bio.SeqUtils import GC, GC123, GC_skew, xGC_skew
from Bio.SeqUtils import MeltingTemp as mt
import numpy as np

GC(seq_3) # GC content 
100-GC(seq_3) # AT content

# Melting point of DNA (in C degrees)
mt.Tm_Wallace(seq_3)
mt.Tm_GC(seq_3)

# GC skew, (G-C)/(G+C)
#(Abundance of GC content in particular regions of the sequence)
# GC skew pos = leading strand
# GC skew neg = lagging strand

#seq_4 = 'AAAATATATAAAAAAAGGGGGGCCCCCCAAAAAAATATATATAAG'
seq_4= 'ATAGACTACTACTAAATATATGTAGTAGTAGTAAATATAGTAAAAAAAGAATGATGGGTGTGTCCTGTCCCCTGAAGTATATATATATATAAAAAAAA'
GC123(seq_4) # Returns total%, and first 3 codons

window=5   # nt span, number of regions will be sequence length / window size rounded up 

skew_vals = GC_skew(seq_4, window=window) #returns GC_skew value within each window (# nt)
plt.plot(list(range(1,len(seq_4), window)),skew_vals)
plt.xlabel('position of first nucleotide within each window (#nt={w})'.format(w=window))
plt.xticks(np.arange(1, len(seq_4), window))
plt.ylabel('(G-C)/(G+C)')
plt.pause(0.5)
plt.show(block=True)

# Subsequences 

from Bio.SeqUtils import nt_search

subseq= 'TC'
print('\n\n search for all subsequence {sub} positions in sequence {seq}: \n'.format(sub=subseq, seq=seq_3))
print(nt_search(str(seq_3), str(subseq)))  
# nt_search needs string input, will not take Seq() in. 

#%% Sequence Alignment
# Search for matches, mismatches and gaps between two sequences

# 1. Global alignment (Needle): best concordance/agreement between all characters 

# 2. Local alignment (Water): focus on subsequences 
# used when query and target sequences differ in length

# BASIC LOCAL ALIGNMENT SEARCH TOOL (BLAST)
 
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

print('\n\n Sequences: \n')
seq_1= Seq('ACTCGT')
seq_2= Seq('ATTCG')

print(seq_1, '\n')
print(seq_2)

# Global alignment
glob_alignments = pairwise2.align.globalxx(seq_1, seq_2)
print('\n\n Possible global alignments: \n')
for a in glob_alignments:
    print(format_alignment(*a))
glob_alignment2 = pairwise2.align.globalxx(seq_1, seq_2, one_alignment_only=True, score_only=True)
    

# Local alignment
loc_alignments = pairwise2.align.localxx(seq_1,seq_2)
print('\n\n Possible local alignments: \n')
for a in loc_alignments:
    print(format_alignment(*a))
loc_alignment2 = pairwise2.align.localxx(seq_1, seq_2, one_alignment_only=True, score_only=True)
    

# Check for similarity (%) using Alignment
# Fraction of nts that are the same / total number of nt (use the larger
# of both sequences)  * 100 
print('\n\nSimilarity between {seq_1} and {seq_2}: \n'.format(seq_1=seq_1, seq_2=seq_2))
print('global alignemnt: ', '{:.2f}'.format(glob_alignment2/len(seq_1)*100), '%') 
print('local alignment: ', '{:.2f}'.format(loc_alignment2/len(seq_1)*100), '%') 

# Similarity is defined as the resemblance between two sequences,
# and comunicates the minimal number of edit operations (indels, etc.) needed
# in order to transform one sequence into an exact copy of the other. 


# Find out all the possible global alignments with the maximum similarity score 
# with parameters: match = +2, mismatch = -1, gap= -0.5, gap_extension= -0.1

print('\n\nFind all the possible global alignments with the maximum similarity score: \n')
glb_alignment = pairwise2.align.globalms(seq_1, seq_2, 2, -1, -0.5,-0.1)
for a in glb_alignment:
    print(format_alignment(*a))
    
# Alignment similarity vs Alignment identity 

# Identity implies the number of characters that match exactly between
# two sequences, not counting gaps. 
# The measure is relational to the shorter sequence, and is not transitive.
# It doesn't translates to both sequences being the exact same either. 


seqA = Seq('AAGGCTT')
seqB = Seq('AAGGC')
seqC = Seq('AAGGCATTT')

AvB = pairwise2.align.localxx(seqA,seqB,one_alignment_only=True, score_only=True)
BvC = pairwise2.align.localxx(seqC,seqB,one_alignment_only=True, score_only=True)
AvC = pairwise2.align.localxx(seqA,seqC,one_alignment_only=True, score_only=True)

print('\n\nIdentity vs Similarity for sequences A ({A}), B ({B}) and C ({C}): \n'.format(A=seqA, B=seqB, C=seqC))

print("AvB", "{:.2f}".format(AvB/len(seqB)*100), "% identical,", "{:.2f}".format(AvB/len(seqA)*100), "% similar ({s} edits needed for exact copy)".format(s=int(len(seqA)-AvB)))
print("BvC", "{:.2f}".format(BvC/len(seqB)*100), "% identical,", "{:.2f}".format(BvC/len(seqC)*100), "% similar ({s} edits needed for exact copy)".format(s=int(len(seqC)-BvC)))
print("AvC", "{:.2f}".format(AvC/len(seqA)*100), "% identical,", "{:.2f}".format(AvC/len(seqC)*100), "% similar ({s} edits needed for exact copy)".format(s=int(len(seqC)-AvC)))


#%% Calculate distance between sequences

seq_1 = Seq('ACTAT')
seq_2 = Seq('ACTTA')
seq_3 = Seq('ACTT')
print('\n\nseq_1=', seq_1, '\nseq_2=', seq_2, '\nseq_1r=', seq_1[::-1], '\nseq_3=', seq_3)

# 1. Hamming Distance: it is done between same-length sequences
# Returns the number of mismatches 

def hamming_distance(lhs,rhs):
    if len(lhs) == len(rhs):
        return len([(x,y) for x,y, in zip (lhs,rhs) if x != y])
    else:
        return print('Hamming Distance Error: unequal lengths')
    
h1v1 = hamming_distance(seq_1,seq_1)
h1v2 = hamming_distance(seq_1,seq_2)
h1v1r = hamming_distance(seq_1, seq_1[::-1])


print('\n\nHamming distance for the same sequence:', h1v1, 'mismatches')
print('\nHamming distance between seq_1 & seq_2:', h1v2, 'mismatches')
print('\nHamming distance between seq_1 & reverse seq_1:', h1v1r, 'mismatches')


# 2. Levenshtein Distance: can be done between unequal length sequences
#  returns the number of edits needed to transform one string into another
# pip install python-Levensthein

from Levenshtein import distance

print('\n\nLevenshtein distance between seq_1 & seq_2=', distance(seq_1, seq_2), 'edition event(s) needed') 
print('\nLevenshtein distance between seq_1 & seq_3=', distance(seq_1, seq_3), 'edition event(s) needed') 


#%% DotPlots for similarity
# helps identify direct or inverted repeats, 
# regions with low complexity
# similar regions
# repeats
# rearrangements
# gene order
# RNA structures
# etc.

seq_1 = Seq('ACTTAG')
seq_2 = Seq('AC')

def delta(x,y):
    return 0 if x == y else 1

def M(seq_1, seq_2, i,j,k):
    return sum(delta(x,y) for x,y in zip(seq_1[i:i+k], seq_2[j:j+k]))

def makeMatrix(seq_1,seq_2,k):
    n = len(seq_1)
    m = len(seq_2)
    return [[M(seq_1,seq_2, i, j, k) for j in range(m-k+1)] for i in range(n-k+1)]

def plotMatrix(M,t,seq_1,seq_2,nonblank=chr(0x25A0), blank=' '):
    print(' |'  + seq_2)
    print('-'*(2+len(seq_2)))
    for label,row in zip(seq_1,M):
        line = ''.join(nonblank if s<t else blank for s in row)
        print(label+'|'+line)

def dotplot(seq_1,seq_2, k=1, t=1):
    M = makeMatrix(seq_1,seq_2, k)
    plotMatrix(M, t, seq_1, seq_2)

# print in console
dotplot(seq_1,seq_2)

# interactive plot using plt
def dotplotx(seq_1,seq_2):
    plt.imshow(np.array(makeMatrix(seq_1,seq_2,1)))
    xt=plt.xticks(np.arange(len(list(seq_2))), list(seq_2))
    yt=plt.yticks(np.arange(len(list(seq_1))), list(seq_1))
    plt.pause(0.5)
    plt.show(block=True)

dotplotx(seq_1,seq_2)

# 

dna1= Seq('ATGATCTCGTAA')
dna2= Seq('ATTATGTCGTAA')

dotplotx(dna1,dna2)

score = pairwise2.align.localxx(dna1,dna2,one_alignment_only=True, score_only=True)

score/len(dna2)*100


#%% FASTA, GENBANK, PDB

# Load FASTA
from Bio import SeqIO

for record in SeqIO.parse("sequence.fasta", "fasta"):
    print(record.id)
    print(record.description)
    
fasta_dna_record = SeqIO.read("sequence.fasta", "fasta")
fasta_dna_seq    = fasta_dna_record.seq 

# Load GenBank
for record in SeqIO.parse("sequence.gb", "genbank"):
    print(record)

gb_dna_record = SeqIO.read("sequence.gb", "gb")
gb_dna_seq    = gb_dna_record.seq 

# do stuff
window=250   # nt span, number of regions will be sequence length / window size rounded up 
skew_vals = GC_skew(gb_dna_seq, window=window) #returns GC_skew value within each window (# nt)
plt.plot(list(range(1,len(gb_dna_seq), window)),skew_vals)
plt.xlabel('position of first nucleotide within each window (#nt={w})'.format(w=window))
plt.xticks(np.arange(1, len(seq_4), window*16))
plt.ylabel('(G-C)/(G+C)')
plt.pause(0.5)
plt.show(block=True)


#%% BLAST

from Bio.Blast import NCBIWWW

# BLAST search for a given sequence using NCBI's QBLAST

# with NCBIWWW.qblast("blastn","nt",fasta_dna_seq[:100]) as result_handle:
#     with open("result_blast_covid.xml","w") as xml_file:
#         xml_file.write(result_handle)

# with open("result_covd_blast.xml", "w") as f:
#     for i in result_handle:
#         f.write(i)
        
    

plt.ioff() 