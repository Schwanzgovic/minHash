
# function to read the genome files and return a list of strings 
# each string should represent one genome 

def readFileAndReturnStringList(filename):
    listOfStrings = []
    string = ""
    # TODO
    return listOfStrings

# this function takes a string and returns all subsequences of this string with length k

def kMers(string, k): 
    listOfKmers = []
    # TODO
    return listOfKmers

# function to compute the min hash for a given hash function and given set of k_mers

def minHashValue(kMers, randomHashfunction):
    # initialize old value with the hash value of the first k-mer
    minHash = randomHashfunction(kMers[0])
    for kMer in kMers:
        newValue = randomHashfunction(kMer)
        if minHash > newValue:
            continue
        else: 
            minHash = newValue

    return minHash

