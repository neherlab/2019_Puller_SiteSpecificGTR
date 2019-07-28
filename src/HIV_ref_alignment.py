import numpy as np
from Bio import AlignIO
from matplotlib import pyplot as plt

def get_distance(aln1, aln2):
    diversity = np.ones(aln1.shape[1])
    for nuc in ['A', 'C', 'G', 'T', '-']:
        diversity -= np.mean(aln1==nuc, axis=0)*np.mean(aln2==nuc, axis=0)
    return diversity

aln = AlignIO.read('HIV_data/HIV1_REF_2010_pol_DNA.fasta', 'fasta')
seq_names = [s.id for s in aln]

aln_array = np.array(aln)

ungapped = np.mean(aln_array=='-', axis=0)<0.2

aln_array_gap_stripped = aln_array[:,ungapped]
overall_div = get_distance(aln_array_gap_stripped, aln_array_gap_stripped)
print("overall diversity:", overall_div.mean())


HIV1_subtypes = ['B', 'C', 'H', 'D']
CPZ_sequences = ['CPZ']

for st1 in HIV1_subtypes:
    ind1 = np.array([s.split('.')[1]==st1 for s in seq_names])
    for st2 in HIV1_subtypes:
        ind2 = np.array([s.split('.')[1]==st2 for s in seq_names])
        dist = get_distance(aln_array_gap_stripped[ind1], aln_array_gap_stripped[ind2])
        print(f"{st1}-{st2}:{dist.mean()}")


for st1 in HIV1_subtypes+CPZ_sequences:
    ind1 = np.array([s.split('.')[1]==st1 for s in seq_names])
    for st2 in CPZ_sequences:
        ind2 = np.array([s.split('.')[1]==st2 for s in seq_names])
        dist = get_distance(aln_array_gap_stripped[ind1], aln_array_gap_stripped[ind2])
        print(f"{st1}-{st2}:{dist.mean()}")





aln_SIV = AlignIO.read('HIV_data/SIV_COM_2017_pol_DNA.fasta', 'fasta')
seq_names_SIV = [s.id for s in aln]

aln_array_SIV = np.array(aln_SIV)
ungapped_SIV = np.mean(aln_array_SIV=='-', axis=0)<0.2
aln_array_gap_stripped_SIV = aln_array_SIV[:,ungapped_SIV]

overall_div_SIV = get_distance(aln_array_gap_stripped_SIV, aln_array_gap_stripped_SIV)
print("overall SIV diversity:", overall_div_SIV.mean())
plt.plot(sorted(overall_div_SIV), np.linspace(0,1,len(overall_div_SIV)))
plt.plot(sorted(overall_div), np.linspace(0,1,len(overall_div)))

