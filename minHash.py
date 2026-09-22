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

# this reads a chosen fasta file, we add it to our final function
#genomes = {rec.id: str(rec.seq) for rec in SeqIO.parse(file_fasta, "fasta")}



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
    m = len(m) 
    #start the sketch with m big start values
    sketch = [float("inf")]*m
    #stream kmers directly from the sequence
    for i in range(len(sequence)-k+1):
        kmer = sequence[i:i+k]
 
        #eveluate kmer against all seeds (streaming)
        for seed in range(m): # we always use the same seeds (0,1,2,3,....,m-1)
            h_val = hash_kmer(kmer, seed)
            if h_val < sketch[seed]:
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
# xxx i dont know if we want to do the gonomic distance as a seperate function or combine it with the jaccard


# this reads a chosen fasta file, we add it to our final function
# genomes = {rec.id: str(rec.seq) for rec in SeqIO.parse(file_fasta, "fasta")}


















