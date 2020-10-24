# PatternCount from Ch1
def PatternCount(Text, Pattern):
    count = 0
    last_position = len(Text) - len(Pattern) + 1

    for x in range(last_position):
        if Text[x : x + len(Pattern)] == Pattern:
            count += 1
    return(count)

#symbol array: counts of e.g. A in a sliding window     # ! slow algorithm won't work for big datasets
def SymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]         # This is needed for circular genomes - n//2 is used becuase we keep track
    for i in range(n):                               #of half the genome's bases from OrI to Ter
        array[i] = PatternCount(ExtendedGenome[i:i+(n//2)], symbol)     # adds no. occurances to value of array{index,value}
    return array

# This prints a dictionary with {key = index, value = no. occrances of A}
print(SymbolArray("AAAAGGGG", 'A'))

#with E.coli genome
Ecoli = open('/Users/panayiotissoteriou/Desktop/panayiotis/online courses/bioinformatics specialisation/Programming meets biology/E_coli_genome.txt', 'r')
Ecoli_genome = Ecoli.read()

print(SymbolArray('Ecoli_genome', 'C'))


# more efficient Symbol array algorithm:
def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]

    # look at the first half of Genome to compute first array value
    array[0] = PatternCount(symbol, Genome[0:n//2])

    for i in range(1, n):
        # start by setting the current array value equal to the previous array value
        array[i] = array[i-1]

        # the current array value can differ from the previous array value by at most 1
        if ExtendedGenome[i-1] == symbol:         # asks if the base that just disappeared out of our moving window was the same as the base we're looking for
            array[i] = array[i]-1                 #if so, we remove one from the number of bases in the current window

        if ExtendedGenome[i+(n//2)-1] == symbol:   # this asks if the base that just came into 'front' of the moving window is the same as the base we're looking for
            array[i] = array[i]+1                  # if so, we add one to the number of bases in the current window
    return array


# Skew array
def SkewArray(Genome):
    Skew = [0]
    for i in range(len(Genome)):
        if Genome[i] == "C":
            Skew.append(Skew[i] - 1)        #skew decreases
        elif Genome[i] == "G":
            Skew.append(Skew[i] + 1)        # skew increases
        else:                               #skew remains unchanged
            Skew.append(Skew[i])

    return Skew

#print(type(SkewArray("CATGGGCATCGGCCATACGCC")))
#print("0 -1 -1 -1 0 1 2 1 1 1 0 1 2 1 0 0 0 0 -1 0 -1 -2")

#minimum Skew
def MinimumSkew(Genome):
    # generate an empty list positions
    positions = list()
    # set a variable equal to SkewArray(Genome)
    array = SkewArray(Genome)
    # find the minimum value of all values in the skew array
    min_array = min(array)
    for i in range(len(array)):
        if min_array == array[i]:
            positions.append(i)
    return positions
    # range over the length of the skew array and add all positions achieving the min to positions


print(MinimumSkew("GATACACTTCCCGAGTAGGTACTG"))

# Hamming distance
def HammingDistance(p, q):
    count = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            count += 1
    return count

#print(HammingDistance("GGGCCGTTGGT", "GGACCGTTGAC"))
print(HammingDistance("CTACAGCAATACGATCATATGCGGATCCGCAGTGGCCGGTAGACACACGT", "CTACCCCGCTGCTCAATGACCGGGACTAAAGAGGCGAAGATTATGGTGTG"))

# approximate Pattern matching
def ApproximatePatternMatching(Text, Pattern, d):
    positions = []
    for i in range(len(Text)-len(Pattern)+1):
        if HammingDistance(Text[i:i+len(Pattern)],Pattern) <= d:
            positions.append(i)
    return positions

print(ApproximatePatternMatching('CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT','ATTCTGGA',3))

#modifying ApproximatePatternMatching to find no. of occurances of k-mers with 1 mismatch
def ApproximatePatternCount(Text, Pattern, d):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if HammingDistance(Text[i:i+len(Pattern)],Pattern) <= d:
            count += 1
    return count

def HammingDistance(p, q):
    count = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            count += 1
    return count

print(ApproximatePatternCount('TTTAGAGCCTTCAGAGG','GAGG', 2))
