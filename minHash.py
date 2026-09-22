import hashlib
from Bio import SeqIO
import zlib

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
genomes = {rec.id: str(rec.seq) for rec in SeqIO.parse(file_fasta, "fasta")}



# this function takes a string and returns all subsequences of this string with length k

#def kMers(string, k): 
 #   listOfKmers = []
  #  for i in range(len(string)-k+1):
   #     listOfKmers.append(string[i:i+k])
   # return listOfKmers

#takes a sequence as string and returns all kmers as intergers in a list (xxx i am using zlib here!)
def get_kmer_ints(sequence,k):
    kmer_ints =[]
    
    #loop through secuenqe and make all kmers
    for i in range(len(sequence) - k +1):
        kmer = sequence[i : i+k]
        
        # converts the string to a unique 32bit interger
        kmer_int = zlib.crc32(kmer.encode('utf-8')) 
        
        kmer_ints.append(kmer_int)
    return kmer_ints

def
# function to get m deterministic hash functions with set seeds
def get_hash_functions(m):
    # TODO
    return #list of hash functions

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
def create_sketch(kMers, hash_functions)
    m = len(hash_functions) 
    #start the sketch with m big start values
    sketch = [float("inf")]*m

    #goes through one k-mer at the time (streaming)
    for kmer in kmers:
        #tests kmer on all m functions
        for i, h in enumerate(hash_functions):
            hash_val = h(kmer)

            #save the value if it is the smallest seen so far
            if hash_val < sketch[i]:
                sketch[i] = hash_val

#takes two lists and compares them element to elemnt
def estimate_jaccard(sketch_A, sketch_B):
    #compares the two sketches, counts how many positions have the same value, divides number of matches with m
    #to give an estimat J(A,B)

    # TODO
    return jaccard_value


# xxx i dont know if we want to do the gonomic distance as a seperate function or combine it with the jaccard

# this reads a chosen fasta file, we add it to our final function
# genomes = {rec.id: str(rec.seq) for rec in SeqIO.parse(file_fasta, "fasta")}







