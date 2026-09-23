from Bio import SeqIO
import zlib
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

# replaces every entry on and below the diagonal with a two
# this code assumes an n x n matrix (square)
# this is to make sure that there is only one unique minimal distance
def lowerTriangleTwos(matrix):
    lowerTriangularTwoMatrix = matrix.copy()
    n = matrix.shape[0]
    for i in range(n):
        for j in range(i+1):
           lowerTriangularTwoMatrix[i][j] = 2 
    return lowerTriangularTwoMatrix

# finds the indices of a matrix where the matrix has a minimal value
def findMinimumIndicesOfMatrix(matrix):
    n = matrix.shape[0]

    # initialize column and rowIndex and the ol minimum with the minimum from the first row
    oldMin = matrix[0][np.argmin(matrix[0])]
    columnIndex = np.argmin(matrix[0])
    rowIndex = 0

    # loop throug the rows of the matrix
    for i in range(n):
        # assign minimal value of that row and safe the its index
        newMinIndex = np.argmin(matrix[i])
        newMin = matrix[i][newMinIndex]
        if newMin < oldMin:
            # if the minimum of the current row is smaller then the old one
            # update row and column index
            columnIndex = newMinIndex
            rowIndex = i
    return rowIndex, columnIndex


# this function takes a distance matrix and a row and colum index 
# which correspnd to two indices of genomes in the genomeID list
# we pick the lower one of these indices and replace it by the common parent of the two
# the higher one gets deleted (so both its distance row and column!)
# (i, j) is the pair with the minimum distance (with i < j)
# identify all remaining indices k (excluding i and j)
def update_matrix(matrix, minIndexA, minIndexB):

    # making sure i < j:
    i = min(minIndexA, minIndexB)
    j = max(minIndexA, minIndexB)
    remaining = [k for k in range(len(matrix)) if k != i and k != j]
    print(remaining)
    # compute the new row of distances from merged node (i,j) to each remaining node k
    new_row = [(matrix[i][k] + matrix[j][k]) / 2.0 for k in remaining]

# 4. Rebuild the distance matrix
    new_matrix = []
    for r_idx, r in enumerate(remaining):
    # Keep existing distances between remaining nodes
        row = [matrix[r][c] for c in remaining]
        # Add distance to the new merged node
        row.append(new_row[r_idx])
        new_matrix.append(row)
    # Add the final row for the new merged node itself (distance to self = 2.0)
    new_matrix.append(new_row + [2.0])
    return np.array(new_matrix)

 
# function to build the evolutionary tree. 
# It returns an array of the form [[parent1, child11, child12], [parent2, child21, child22], ...]
# from this list, the tree can be constructed bottom-up by going through the list from left to right and 
# an assigning the parent to the two children
def constructEvolutionaryTree(distance_matrix, genome_IDS):
    n = len(genome_IDS)
    # initilize updateMatrixForIndexSearch and the updateMatrix (which is the updated distance matrix)
    updatedMatrixForIndexSearch = lowerTriangleTwos(distance_matrix)
    updateMatrix = distance_matrix.copy()
    # initialize the evolution tree as empty list 
    evolutionaryTree = []
    # copy the genome IDs and make sure the mainupalion of the does not affect the original IDs
    treeIDs = genome_IDS.copy()
    # the loop has to run exactly n-1 times
    # since exactly n-1 parents can be found, where the n-1th is the root 
    for k in range(n-1): 
        # extract the indices corresponding to the minimum distance for an unspecified genome A and genome B 
        minIndexA, minIndexB = findMinimumIndicesOfMatrix(updatedMatrixForIndexSearch)
        # update the distance matrixwith these indices
        updateMatrix = update_matrix(updateMatrix, minIndexA, minIndexB)
        # take the now updated distance matrix and make it ready for the next index search
        updatedMatrixForIndexSearch = lowerTriangleTwos(updateMatrix)

        # assign the smaller index to i
        i = min(minIndexA, minIndexB)
        # and the bigger index to j
        j = max(minIndexA, minIndexB)

        # the "parent" nodes (or genomes) are identified by increasing ascii symbolds starting from 'A'
        parentNode = chr(65 + k)

        # extract the childs from the tree ID list (note that we do not use the genome ID list here
        # since it also needs to contain the added parents)
        childOne = treeIDs[i]
        childTwo = treeIDs[j]

        # add the [parent, child1, child2] triple to the evolution tree
        evolutionaryTree.append([parentNode, childOne, childTwo])
        # update the tree IDs
        treeIDs[i] = parentNode
        # remove the bigger index from the list 
        treeIDs.remove(childTwo)
    return evolutionaryTree
    
#Solves Task 1    
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


#Solves Task 2
#creates the evolution tree from the function
evolutionaryTree = constructEvolutionaryTree(dist_matrix, genome_ids)


# plots the list, which remember can be built bottom up
print("evolution tree: ", "\n", evolutionaryTree)





