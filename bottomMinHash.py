from Bio import SeqIO
import zlib
import hashlib
import random
import numpy as np


def hash_kmer(kmer, seed):
  # crc32 takes a seed and a kmer and hashes 
  return zlib.crc32(kmer.encode("utf-8"), seed)


    
def create_bottom_sketch(sequence, k, m):
    #creates a bottom hash sketch by only running one hashfunction
    #and saving the m smallest hash values
    
    seen_hashes = set()
    
    for i in range(len(sequence)-k+1):
        kmer = sequence[i : i+k]
        h_val = hash_kmer(kmer, seed = 2001)
        seen_hashes.add(h_val)
    
    return sorted(seen_hashes)[:m]


def estimate_jaccard_bottom(sketch_A, sketch_B):
    m=len(sketch_A)
    
    set_A = set(sketch_A)
    set_B = set(sketch_B)
    
    #join the sets and take the ma smallest
    union_m = sorted(set_A.union(set_B))[:m]
    
    #count how many of the m smallest are in both A and B
    matches = sum(1 for h in union_m if h in set_A and h in set_B)
    
    return float(matches/m)


#takes a dict of sketches
#calculates a NxN distance (1-Jaccard) as a numpy-array
def create_distance_martix(sketches):
    
    genome_ids = list(sketches.keys())
    N = len(genome_ids)
    
    matrix = np.zeros((N,N))
    
    #calc all pairwise distances
    for i in range(N):
        for j in range(i+1,N):
            sim = estimate_jaccard_bottom(sketches[genome_ids[i]], sketches[genome_ids[j]])
            dist = 1 - sim
            
            #symmetric
            matrix[i,j] = dist
            matrix[j,i] = dist
    
    return matrix, genome_ids
    
#Solves problem 1    
# 1. read all genomes from FASTA-file
genomes = {rec.id: str(rec.seq) for rec in SeqIO.parse("ecoli_20_genomes_short_acc.fa", "fasta")} #change "test3.fa" to wanted file

# 2. choose parameters
k = 21  # standard for bacteria (i googled, but we might want to play around with it and find a justification)
m = 100  # number of hash functions / seeds

# 3. create sketches-dictionary
sketches = {}


for genome_id, sequence in genomes.items():
  # calls create_sketch-function for every genom
  sketches[genome_id] = create_bottom_sketch(sequence, k, m)

print(f"Created {len(sketches)} sketches") # xxx just to check we can remove later.

#print(sketches)
dist_matrix, genome_ids = create_distance_martix(sketches) 

print("Distance matrix: ", "\n", dist_matrix,"\n", "Genomde IDs: ", "\n", genome_ids)  
    

