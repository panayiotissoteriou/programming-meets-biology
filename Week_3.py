print('Hello world!')

#the frequent words algorothm can identify most frequent k-mers
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

                    #remove Enters
print(FrequentWords('atgaccgggatactgataaaaaaaagggggggggcgtacacattagataaacgtatgaagtacgttagactcggcgccgccgacccctattttttgagcag
atttagtgacctggaaaaaaaatttgagtacaaaacttttccgaataaaaaaaaagggggggatgagtatccctgggatgacttaaaaaaaagggggggtgctctcccgatttttgaatatg
taggatcattcgccagggtccgagctgagaattggatgaaaaaaaagggggggtccacgcaatcgcgaaccaacgcggacccaaaggcaagaccgataaaggagatcccttttgcggtaatg
tgccgggaggctggttacgtagggaagccctaacggacttaataaaaaaaagggggggcttataggtcaatcatgttcttgtgaatggatttaaaaaaaaggggggggaccgcttggcgcac
ccaaattcagtgtgggcgagcgcaacggttttggcccttgttagaggcccccgtaaaaaaaagggggggcaattatgagagagctaatctatcgcgtgcgtgttcataacttgagttaaaaa
aaagggggggctggggcacatacaagaggagtcttccttatcagttaatgctgtatgacactatgtattggcccattggctaaaagcccaacttgacaaatggaagatagaatccttgcata
aaaaaaagggggggaccgaaagggaagctggtgagcaacgacagattcttacgtgcattagctcgcttccggggatctaatagcacgaagcttaaaaaaaaggggggga', 15))


# Count(motifs):
def Count(Motifs):
    count = {}
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(0)

    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count

print(Count(['AACGTA','CCCGTT','CACCTT','GGATTA','TTCCGG']))

#profile (motifs)

def Profile(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}
    # insert your code here
    profile1 = Count(Motifs)
    for key in "ACGT":
        for key in profile1:
            profile[key] = [float(x) / float(t) for x in profile1[key]]
    return profile


print(Profile(['AACGTA','CCCGTT','CACCTT','GGATTA','TTCCGG']))

# consensus(motif)
def Count(Motifs):
    count = {}
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(0)

    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count

def Consensus(Motifs):
    k = len(Motifs[0])
    count = Count(Motifs)
    consensus = ""

    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol

    return consensus

print(Consensus(['AACGTA','CCCGTT','CACCTT','GGATTA','TTCCGG']))

# Score(motifs)
def Consensus(Motifs):
    k = len(Motifs[0])
    count = Count(Motifs)
    consensus = ""

    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol

    return consensus
def Count(Motifs):
    count = {}
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(0)

    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count

def Score(Motifs):
    consensus = Consensus(Motifs)
    count = 0
    for motif in Motifs:
        for i in range(len(motif)):
            if motif[i] != consensus[i]:
                count += 1
    return count

print(Score(['AACGTA','CCCGTT','CACCTT','GGATTA','TTCCGG']))

#Pr function
def Pr(Text, Profile):
    p = 1
    for i in range(len(Text)):
        p = p * Profile[Text[i]][i]
    return p

Profile={"A": [0.4, 0.3, 0.0, 0.1, 0.0, 0.9],
"C": [0.2, 0.3, 0.0, 0.4, 0.0, 0.1],
"G":[0.1, 0.3, 1.0, 0.1, 0.5, 0.0],
"T":[0.3, 0.1, 0.0, 0.4, 0.5, 0.0]}

print(Pr("AAGTTC", Profile))


print(Pr('ACGGGGATTACC', [0.2 0.2 0.0 0.0 0.0 0.0 0.9 0.1 0.1 0.1 0.3 0.0
0.1 0.6 0.0 0.0 0.0 0.0 0.0 0.4 0.1 0.2 0.4 0.6
0.0 0.0 1.0 1.0 0.9 0.9 0.1 0.0 0.0 0.0 0.0 0.0
0.7 0.2 0.0 0.0 0.1 0.1 0.0 0.5 0.8 0.7 0.3 0.4]))

# profile most probable pattern
def Pr(Text, Profile):
    p = 1
    for i in range(len(Text)):
        p = p * Profile[Text[i]][i]
    return p

def ProfileMostProbablePattern(text, k, profile):
    p=-1
    result=text[0:k]
    for i in range(len(text)-k+1):
        seq=text[i:i+k]
        pr=Pr(seq,profile)
        if pr>p:
            p=pr
            result=seq
    return result

#GreedyMotifSearch
def Count(Motifs):
    count = {}
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(0)

    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count

def Consensus(Motifs):
    k = len(Motifs[0])
    count = Count(Motifs)
    consensus = ""

    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol

    return consensus

def Profile(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}
    # insert your code here
    profile1 = Count(Motifs)
    for key in "ACGT":
        for key in profile1:
            profile[key] = [float(x) / float(t) for x in profile1[key]]
    return profile

def Score(Motifs):
    consensus = Consensus(Motifs)
    count = 0
    for motif in Motifs:
        for i in range(len(motif)):
            if motif[i] != consensus[i]:
                count += 1
    return count

def Pr(Text, Profile):
    p = 1
    for i in range(len(Text)):
        p = p * Profile[Text[i]][i]
    return p

def ProfileMostProbablePattern(text, k, profile):
    p=-1
    result=text[0:k]
    for i in range(len(text)-k+1):
        seq=text[i:i+k]
        pr=Pr(seq,profile)
        if pr>p:
            p=pr
            result=seq
    return result

def GreedyMotifSearch(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])

    n = len(Dna[0])

    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs

    return BestMotifs
