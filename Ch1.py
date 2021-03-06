a = "sup"
print(a)

Vibrio_chol_genome = 'ATCAATGATCAACGTAAGCTTCTAAGCATGATCAAGGTGCTCACACAGTTTATCCACAACCTGAGTGGATGACATCAAGATAGGTCGTTGTATCTCCTTCCTCTCGTACTCTCATGACCACGGAAAGATGATCAAGAGAGGATGATTTCTTGGCCATATCGCAATGAATACTTGTGACTTGTGCTTCCAATTGACATCTTCAGCGCCATATTGCGCTGGCCAAGGTGACGGAGCGGGATTACGAAAGCATGATCATGGCTGTTGTTCTGTTTATCTTGTTTTGACTGAGACTTGTTAGGATAGACGGTTTTTCATCACTGACTAGCCAAAGCCTTACTCTGCCTGACATCGACCGTAAATTGATAATGAATTTACATGCTTCCGCGACGATTTACCTCTTGATCATCGATCCGATTGAAGATCTTCAATTGTTAATTCTCTTGCCTCGACTCATAGCCATGATGAGCTCTTGATCATGTTTCCTTAACCCTCTATTTTTTACGGAAGAATGATCAAGCTGCTGCTCTTGATCATCGTTTC'

print(Vibrio_chol_genome)
print(len(Vibrio_chol_genome))          #540
print(range(len(Vibrio_chol_genome)))   #[0, 1, 2, ..., 539]: shows indexing

#counting patterns
#sliding window:
    # for i in range(len(Text)-len(Pattern)+1):

kmer_pattern = 'ATCAATGATC'  # k-mer length = 10

count = 0
last_position = len(Vibrio_chol_genome) - len(kmer_pattern) +1
for i in range(last_position):                          # for every element between 0 and 531 (540-10+1) [+1 bc 0-indexing is used]
    if Vibrio_chol_genome[i : i + len(kmer_pattern)] == kmer_pattern:     # analyse if elements in positions 0 to 530 (540-len(kmer)) same as pattern
        count = count + 1
print(count)

#Pattern count function
def patterncount(Text, Pattern):
    count = 0
    last_position = len(Text) - len(Pattern) + 1

    for x in range(last_position):
        if Text[x : x + len(Pattern)] == Pattern:
            count += 1
    return(count)

print(patterncount(Vibrio_chol_genome, 'ATC'))

# Patterncount 2: (one to use in course)
def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return(count)

# exercise: how many times 'TGATCA' occurs in the vibrio genome
PatternCount(Vibrio_chol_genome, 'TGATCA')

# frequency map function:


def FrequencyMap(Text, k):
    freq = {}                   #emtpy dict
    n = len(Text)
    # adds pattern into the dictionary
    for i in range(n-k+1):      #at every index from 0 to the last pattern of k-length
        Pattern = Text[i:i+k]   # pattern = sliding kmer window
        freq[Pattern] = 0       # adds pattern into the dict

    # counts pattern occurances
        for i in range(n-k+1):          # at every index from 0 to the last pattern of k-length
            if Text[i:i+k] == Pattern:  # if sequences at all possible indexes are the same as the pattern
                freq[Pattern] = freq.get(Pattern, 0) + 1        # add 1 to the dictionary value corresponding to the key
    return freq

# example:
print(FrequencyMap('ACGTGTGTGTGTGGACGATGACTG', 3))

# list of most frequent words

def FrequencyMap(Text, k):
    freq = {}                   #emtpy dict
    n = len(Text)
    # adds pattern into the dictionary
    for i in range(n-k+1):      #at every index from 0 to the last pattern of k-length
        Pattern = Text[i:i+k]   # pattern = sliding kmer window
        freq[Pattern] = 0       # adds pattern into the dict

    # counts pattern occurances
        for i in range(n-k+1):          # at every index from 0 to the last pattern of k-length
            if Text[i:i+k] == Pattern:  # if sequences at all possible indexes are the same as the pattern
                freq[Pattern] = freq.get(Pattern, 0) + 1        # add 1 to the dictionary value corresponding to the key
    return freq

def FrequentWords(Text, k):
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        if freq[key] == m:
            pattern = key
            words.append(pattern)
    return words

#example
#print(FrequentWords("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4))
print(FrequentWords("TAAACGTGAGAGAAACGTGCTGATTACACTTGTTCGTGTGGTAT", 3))

# Reverse complementary strand
text = "ACCCGTTGTGTGACG"
comp = {"A":"T", "T":"A", "C":"G", "G":"C"}
complement = ""
for i in text:
    complement += comp[i]
    print(complement[::-1])

def ReverseComplement(Text):
    comp = {"A":"T", "T":"A", "C":"G", "G":"C"}
    complement = ""
    for i in Text:
        complement += comp[i]
    return complement[::-1]

print(ReverseComplement("ACGTGGTGTGA"))

# Reverse() function
# using slicing
str = "GTGCGAAC"
stringLen = len(str)
slicedStr = str[stringLen::-1]
print(slicedStr)

#using str.join() & str.reversed()
text = "GTGCGAAC"
reversed_text = "".join(reversed(text))
print(reversed_text)

# the shortest way:
Pattern = "ACGTTTG"
reversed_Pattern=Pattern[::-1]
print(reversed_Pattern)


# prints the complement
text = "GAGACT"
comp = {"G":"C", "C":"G", "T":"A", "A":"T"}
complement_text = ""
for i in text:
    complement_text += comp[i]
print(complement_text)

#to print reverse complement:
text = "GAGACT"
comp = {"G":"C", "C":"G", "T":"A", "A":"T"}
complement_text = ""
for i in text:
    complement_text += comp[i]
print(complement_text[::-1])        #Here is only difference

# to print reverse:
text = "GAGACT"
reverse = ""
for i in text:
    reverse += i
print(reverse[::-1])



# Pattern matching
def PatternMatching(Pattern, Genome):
    positions = []
    last_position = len(Genome) - len(Pattern) + 1

    for i in range(last_position):
        if Pattern == Genome[i:i+len(Pattern)] :
            positions.append(i)
    return positions

print(PatternMatching('ATAT', 'GATATATGCATATACTT'))

