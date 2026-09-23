import hashlib
from Bio import SeqIO
import zlib
import random
import numpy as np


def hash_kmer(kmer, seed):
  # combines seed and kmer into a deterministic 32-bit integer (0 to 2^32 - 1)
  # N: but it is pseudo random?
  hash_input = f"{seed}:{kmer}".encode("utf-8")
  # is this pseudo random?
  return int(hashlib.md5(hash_input).hexdigest()[:8], 16)



# this function creates a sketch (xxx like shown in the lecture : streaming) from a list of kmers and a list of m different hash functions
def create_sketch(sequence, k, m):
    #start the sketch with m big start values
    # N: is it guaranteed here that our hash-values will be lower then this big number??
    # N: since we take 
    sketch = [float("inf")]*m
    #stream kmers directly from the sequence
    for i in range(len(sequence)-k+1):
        kmer = sequence[i:i+k]
 
        #eveluate kmer against all seeds (streaming)
        for seed in range(m): # we always use the same seeds (0,1,2,3,....,m-1)
            h_val = hash_kmer(kmer, seed) # N: note that the seeds here are a hash function
            if h_val < sketch[seed]: # N: compare with the hash value from previous k-mer with same hash-function 
 #MinHash
                sketch[seed] = h_val
    return sketch

 
#takes two lists and compares them element to elemnt
def estimate_jaccard(sketch_A, sketch_B):
	#compares the two sketches, counts how many positions have the same value, divides number of matches with m
	#to give an estimate J(A,B)

	# TODO
	#print(sketch_A)
	#print(sketch_B)
	match = 0
	for i in range(len(sketch_A)):
		if(sketch_A[i]==sketch_B[i]):
			match=match+1
	return float(match/len(sketch_A))


#takes a dict of sketches
#calculates a NxN distance (1-Jaccard) as a numpy-array
def create_distance_martix(sketches):
    genome_ids = list(sketches.keys())
    N = len(genome_ids)
    
    matrix = np.zeros((N,N))
    
    #calc all pairwise distances
    for i in range(N):
        for j in range(i+1,N):
            sim = estimate_jaccard(sketches[genome_ids[i]], sketches[genome_ids[j]])
            dist = 1 - sim
            
            # since matrix is symmetric
            matrix[i,j] = dist
            matrix[j,i] = dist
    return matrix, genome_ids

# find the root of the tree 
# the idea is to take the avergae distance to all other nodes. The root should have the shortest one
# this function returns the id of the genome which most likely is the root 
def findRoot(averageDistances, genomeIDs):
     # get minimal index 
     minIndex = np.argmin(averageDistances)
     return genomeIDs[minIndex]


# function to compute the average distances between each genome and all others
def computeAverageDistances(dist_matrix):
    n = len(dist_matrix[0][:])
    averageDistances = np.zeros(n)
    for i in range(n): 
        sum = 0
        for j in range(n):
            if i!=j:
                sum += dist_matrix[i][j]
        averageDistances[i] = sum/(n-1)
    return averageDistances

# replaces every entry on and below the diagonal with a two
# this code assumes an n x n matrix (square)
def lowerTriangleTwos(matrix):
    lowerTriangularTwoMatrix = matrix
    n = matrix.shape[0]
    for i in range(n):
        for j in range(i+1):
           matrix[i][j] = 2 
    return np.array(lowerTriangularTwoMatrix)

# finds the indices of a matrix where the matrix has a minimal value
def findMinimumIndicesOfMatrix(matrix):
    n = matrix.shape[0]
    oldMin = matrix[0][np.argmin(matrix[0])]
    columnIndex = np.argmin(matrix[0])
    rowIndex = 0
    for i in range(n):
        newMinIndex = np.argmin(matrix[i])
        newMin = matrix[i][newMinIndex]
        if newMin < oldMin:
            columnIndex = newMinIndex
            rowIndex = i
    return rowIndex, columnIndex

# estimates the distances from the parent node to all other nodes by taking the average distances of its children

def estimateParentDistances(distanceMatrix, minIndexA, minIndexB):
    n = distanceMatrix.shape[0]
    parentDistances = np.ones(n-1)
    minAB = min(minIndexA, minIndexB)
    maxAB = max(minIndexA, minIndexB)
    parentDistances[minAB] = 0
    for k in range(n):
        if k == maxAB:
            continue
        elif k != minAB:
            parentDistances[k] = (distanceMatrix[k][minIndexA] + distanceMatrix[minIndexB])/2
        
    return parentDistances        

def updateDistanceMatrix(distanceMatrix, minIndexA, minIndexB):
    n = distanceMatrix.shape[0]
    updateMatrix = np.ones((n-1, n-1))
    minAB = min(minIndexA, minIndexB)
    maxAB = max(minIndexA, minIndexB)
    # TODO
    for i in range(n-1):
        for j in range(n-1):
            if i < minAB and j < minAB:
                updateMatrix[i][j] = distanceMatrix[i][j]

            elif i != minAB and j == minAB:
                # for the first of the two indices, put the parent node (or ancestor here)
                # we estimate it's distances to all other nodes by taking the average of the children
                updateMatrix[i][j] = (distanceMatrix[i][j] + distanceMatrix[i][maxAB])/2
            elif i == minAB and j != minAB:
                updateMatrix[i][j] = 1
            elif i < maxAB and j < maxAB:
                # when both indices are between the minIndexA and minIndexB just copy the entries
                updateMatrix[i][j] = distanceMatrix[i][j]

            elif i < maxAB and j >= maxAB:
                # skip the maxAB column
                updateMatrix[i][j] = distanceMatrix[i][j+1]
            elif i >= maxAB and j < maxAB:
                # skip the maxAB row
                updateMatrix[i][j] = distanceMatrix[i+1][j]
            else:
                # skip the maxAB row and maxAB column
                updateMatrix[i][j] = distanceMatrix[i+1][j+1]

    return updateMatrix

 
# function to build the evolutionary tree. It returns an array which satisfies if paranet is at array[n]
# then array[2n+1] and array[2n+2] are its children
def constructEvolutionaryTree(distance_matrix, genome_IDS):
    n = len(genome_IDS)
    # returns the upper triangle of matrix and cuts out the zeros (also cuts out the doubled distances)
    updatedMatrix = lowerTriangleTwos(distance_matrix)
    evolutionaryTree = []
    # the loop has to run exactly n-2 times
    for k in range(n-1): 
        minIndexA, minIndexB = findMinimumIndicesOfMatrix(updatedMatrix)
        updatedMatrix = updateDistanceMatrix(updatedMatrix, minIndexA, minIndexB)

        parentNode = chr(65 + k)
        childOne = genome_IDS(minIndexA)
        childTwo = genome_IDS(minIndexB)
        evolutionaryTree.append([parentNode, childOne, childTwo])

    return evolutionaryTree
#Solves problem 1    
# 1. read all genomes from FASTA-file
genomes = {rec.id: str(rec.seq) for rec in SeqIO.parse("test1.fa", "fasta")} #change "test3.fa" to wanted file

# 2. choose parameters
k = 21  # standard for bacteria (i googled, but we might want to play around with it and find a justification)
m = 1000  # number of hash functions / seeds (N: can be pretty high for "test1.fa")

# 3. create sketches-dictionary
sketches = {}

for genome_id, sequence in genomes.items():
  # calls create_sketch-function for every genom
  sketches[genome_id] = create_sketch(sequence, k, m)

print(f"Created {len(sketches)} sketches") # xxx just to check we can remove later.

dist_matrix, genome_ids = create_distance_martix(sketches) 

print("Distance matrix: ", "\n", dist_matrix,"\n", "Genomde IDs: ", "\n", genome_ids)
triangularMatrix = lowerTriangleTwos(dist_matrix)
minIndexRow, minIndexColumn = findMinimumIndicesOfMatrix(triangularMatrix)

updateMatrix = updateDistanceMatrix(dist_matrix, minIndexRow, minIndexColumn)
averageDistancesToAllOtherNodes = computeAverageDistances(dist_matrix)
parentDistances = estimateParentDistances(dist_matrix, minIndexRow, minIndexColumn)
#print("average Distances: ", averageDistancesToAllOtherNodes)
print("Upper traingle of matrix: ", "\n", triangularMatrix)
print("Min Index of upper triangle: ", minIndexRow, " ", minIndexColumn)
print("Distances for the parent node: ", "\n", parentDistances)
#print("And finally the updated Matrix: ", "\n", updateMatrix)
