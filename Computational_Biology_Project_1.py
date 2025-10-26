# %%
import collections
from collections import Counter
from Bio import SeqIO
import matplotlib.pyplot as plt
from scipy.stats import binom
import numpy as np

# %%
# Get k-mers and their coverage
def count_fasta_kmers(filename, k):
    counts=collections.defaultdict(int) # Creates empty dict
    N_r=0
    for record in SeqIO.parse(filename, "fasta"): 
        read=str(record.seq)
        N_r+=1
        num_kmers=len(read)-k+1
        for i in range(num_kmers):
            kmer=read[i:i+k]
            counts[kmer]+=1
    return counts, N_r

def get_axes(counts):
    kmer_counts=list(counts.values()) # Extract coverage of each k-mer

    # Counts how many k-mers have a certain coverage
    freq_of_counts=Counter(kmer_counts) # Creates a dictionary where keys=counts.values

    x=sorted(freq_of_counts.keys())
    y=[freq_of_counts[i] for i in x] # List of sorted values

    return x, y

def plot_kmer_frequency_distribution(x,  y, logx=False, logy=False):
    total_kmers = sum(y)
    empirical_prob = np.array(y) / total_kmers

    plt.figure(figsize=(9,6))
    plt.bar(x, empirical_prob, width=0.8, align='center', alpha=0.8) # Creates hist

    plt.xlabel("k-mer coverage")
    plt.ylabel("k-mers with \"x\" coverage")
    plt.title("Distribution of k-mers per coverage")

    if logx:
        plt.xscale('log')
    if logy:
        plt.yscale('log')

    plt.tight_layout()
    plt.show()

def plot_kmer_with_theoretical(x, y, N_r, pX_k, pY_k):
    
    pmfX = binom.pmf(x, N_r, pX_k)
    pmfY = binom.pmf(x, N_r, pY_k)

    plt.figure(figsize=(9,6))
    plt.plot(x, pmfX, 'r-', lw=2, label="Binomial X_k")
    plt.plot(x, pmfY, 'g--', lw=2, label="Binomial Y_k")
    plt.xlabel("k-mer coverage")
    plt.ylabel("Probability")
    plt.title("Theoretical binomial distributions")
    plt.legend()
    plt.tight_layout()
    plt.show()

# %%
# Code with right K

# File paths
fasta_file_wt="salmonella-enterica.reads.fna"
fasta_file_m="salmonella-enterica-variant.reads.fna"

# %%
# Length of the k-mer
k=31

# Sequencing error
nu=0.01

# Fix the genome size
genome_size=5e6
read_length=250

# Probability for a k-mers to overlap a read
p=(read_length-k)/(genome_size-read_length)

# Compute the empirical probabilities for X_k and Y_k
pX_k=p*(1-nu)**k
pY_k =p*nu*(1-nu)**(k-1)

# Build dataset of k-mers and their coverage and take the number of reads
counts_wt, N_r=count_fasta_kmers(fasta_file_wt, k)
x, y=get_axes(counts_wt)

# Hist of empirical distribution
plot_kmer_frequency_distribution(x, y, False, False)

plot_kmer_with_theoretical(x, y, N_r, pX_k, pY_k)

# %%
# Extract solid k-mers

threshold=17.5

solid_kmers_wt=[]
for key, v in counts_wt.items():
    if v > threshold:
        solid_kmers_wt.append(key)

# %%
# Do the same with mutated one

# Build dataset of k-mers and their coverage and take the number of reads
counts_m, _=count_fasta_kmers(fasta_file_m, k)

solid_kmers_m = []
for key, v in counts_m.items():
    if v > threshold:
        solid_kmers_m.append(key)

# %%
# Find k-mers wiith mutation

# Find intersection
intersection = set(solid_kmers_wt) & set(solid_kmers_m)
not_matched_wt = set(solid_kmers_wt) - intersection
not_matched_m  = set(solid_kmers_m)  - intersection

# %%
# Mask based search algorithm in order to find the mutation

def build_mask_index(kmer_list):
    index = {}
    for km in kmer_list:
        for i in range(len(km)):
            mask = km[:i] + '*' + km[i+1:]
            if mask not in index:
                index[mask] = []
            index[mask].append(km)
    return index

def find_one_base_diff(kmers_a, kmers_b):
    index_a = build_mask_index(kmers_a)
    pairs = []
    for km in kmers_b:
        for i in range(len(km)):
            mask = km[:i] + '*' + km[i+1:]
            if mask in index_a:
                for match in index_a[mask]:
                    if sum(aa != bb for aa, bb in zip(km, match)) == 1:
                        pairs.append((match, km))
    return pairs

# %%
# Find and print mutations
pairs = find_one_base_diff(not_matched_wt, not_matched_m)
print(f"{len(pairs)} k-mers differs by one basis.")
print(pairs) #len(pairs) #4     #len(pairs[0]) #2

# %%
# Recover bigger sequence (read)

sequences=[]

for i in range(len(pairs)):
    for j in SeqIO.parse(fasta_file_m, "fasta"):
        read=str(j.seq)
        if pairs[i][1] in read:
            sequences.append(read)
            break
        
print(sequences)