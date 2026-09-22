import hashlib
from Bio import SeqIO
import zlib
import random


# function to read the genome files and return a list of strings (xxx note: maybe we want it as a dictionary, thinking about problem 2)
# each string should represent one genome 

#def readFileAndReturnStringList(filename):
#    listOfStrings = []
#    string = ""
#	with open(filename, 'r') as file:
#		for line in file:
#			if(line[0]=='>' and string != ''):
#				listOfStrings.append(string)
#				string=""
#			else:
#				if(line[0]!='>'):
#					string += line.strip()
#	listOfStrings.append(string)
#	return listOfStrings


# this function takes a string and returns all subsequences of this string with length k
#def kMers(string, k): 
 #   listOfKmers = []
  #  for i in range(len(string)-k+1):
   #     listOfKmers.append(string[i:i+k])
    #return listOfKmers


def hash_kmer(kmer, seed):
  # combines seed and kmer into a deterministic 32-bit integer (0 to 2^32 - 1)
  hash_input = f"{seed}:{kmer}".encode("utf-8")
  return int(hashlib.md5(hash_input).hexdigest()[:8], 16)


# function to compute the min hash for a given hash function and given set of k_mers
#def minHashValue(kMers, randomHashfunction):
    # initialize old value with the hash value of the first k-mer
#    minHash = randomHashfunction(kMers[0])
 #   for kMer in kMers:
  #      newValue = randomHashfunction(kMer)
   #     if minHash > newValue: # xxx i think it should be < but not sure so i keep it as a comment
    #        continue
     #   else: 
      #      minHash = newValue
   # return minHash


# this function creates a sketch (xxx like shown in the lecture : streaming) from a list of kmers and a list of m different hash functions
def create_sketch(sequence, k, m):
    #start the sketch with m big start values
    sketch = [float("inf")]*m
    #stream kmers directly from the sequence
    for i in range(len(sequence)-k+1):
        kmer = sequence[i:i+k]
 
        #eveluate kmer against all seeds (streaming)
        for seed in range(m): # we always use the same seeds (0,1,2,3,....,m-1)
            h_val = hash_kmer(kmer, seed)
            if h_val < sketch[seed]: #MinHash
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
    
#Solves problem 1    
# 1. read all genomes from FASTA-file
genomes = {rec.id: str(rec.seq) for rec in SeqIO.parse("test3.fa", "fasta")} #change "test3.fa" to wanted file

# 2. choose parameters
k = 21  # standard for bacteria (i googled, but we might want to play around with it and find a justification)
m = 100  # number of hash functions / seeds

# 3. create sketches-dictionary
sketches = {}

for genome_id, sequence in genomes.items():
  # calls create_sketch-function for every genom
  sketches[genome_id] = create_sketch(sequence, k, m)

print(f"Created {len(sketches)} sketches") # xxx just to check we can remove later.

dist_matrix, genome_ids = create_distance_martix(sketches) 

print(dist_matrix, genome_ids)   
# this reads a chosen fasta file, we add it to our final function
# genomes = {rec.id: str(rec.seq) for rec in SeqIO.parse(file_fasta, "fasta")}




# this reads a chosen fasta file, we add it to our final function
#genomes = {rec.id: str(rec.seq) for rec in SeqIO.parse(file_fasta, "fasta")}


#PROBLEM: wht k to use
#takes a list of genomes and returns a "good" k value (3 to 600) XXX I have no idea what a good k value is??
def get_good_k(genomes):
	maxlen = 0
	for i in genomes:
		if(len(i)>maxlen):
			maxlen=len(i)
	#print(int(maxlen/14))
	k = int(maxlen/14)
	if(k>600):
		k=600
	#print(k)
	return k

#print(readFileAndReturnStringList('lessontest.txt'))
#gen = readFileAndReturnStringList("ecoli_20_genomes_short_acc.fa")
#gen = readFileAndReturnStringList("lessontest.txt")
gen = readFileAndReturnStringList("test1.fa")

k = get_good_k(gen)
hashesneeded = 0
kmers = []
for i in gen:
	kmers.append(kMers(i,k))
for i in kmers:
	if(len(i)>hashesneeded):
		hashesneeded=len(i)
h = get_hash_functions(hashesneeded)
#print(kmers[0])
#create_sketch(kmers[0],h)
sketches = []
for i in range(len(gen)):
	#print(kmers[i])
	sketches.append(create_sketch(kmers[i],h))


















