print('hello world!')

#countwithpseudocounts()
def CountWithPseudocount(Motifs):
        count = {}
        t = len(Motifs)
        k = len(Motifs[0])
        for symbol in "ACGT":
            count[symbol] = []
            for j in range(k):
                 count[symbol].append(1)


        for i in range(t):
            for j in range(k):
                symbol = Motifs[i][j]
                count[symbol][j] += 1
        return count

#test_data = (\
    'AACGTA',
    'CCCGTT',
    'CACCTT',
    'GGATTA',
    'TTCCGG')

#print(CountWithPseudocount(test_data))

# Profile with Pseudocounts
def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    cont=CountWithPseudocounts(Motifs)
    for i in range(k):
        su=0
        for symbol in "ACGT":
            su=su+cont[symbol][i]
        for symbol in "ACGT":
            cont[symbol][i] = cont[symbol][i]/su
    profile=cont
    return profile

# Input: A set of kmers Motifs
# Output: CountWithPseudocounts(Motifs)
# HINT: You need to use CountWithPseudocounts as a subroutine of ProfileWithPseudocounts
def CountWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    count = {} # output variable
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(1)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count

test_data = (\
    'AACGTA',
    'CCCGTT',
    'CACCTT',
    'GGATTA',
    'TTCCGG')

print(ProfileWithPseudocounts(test_data))

# GreedyMotifSearchWithPseudocounts(Dna, k, t)
def GreedyMotifSearchWithPseudocounts(Dna,k,t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for m in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][m:m+k])
        for j in range(1, t):
            P = ProfileWithPseudocounts(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

def ProfileWithPseudocounts(Motifs):
    profile = {}
    k = len(Motifs[0])
    t = len(Motifs)
    profile = CountWithPseudocounts(Motifs)

    for i in range(k):
        p = 0

        for symbol in "ACGT":
            p = p+profile[symbol][i]
        for symbol in "ACGT":
            profile[symbol][i] = profile[symbol][i]/p
    return profile

def Score(Motifs):
    consensus = Consensus(Motifs)
    count = CountWithPseudocounts(Motifs)

    k = len(Motifs[0])
    t = len(Motifs)
    c = 0
    for Motif in Motifs:
        for i in range(k):
                if Motif[i] != consensus[i]:
                    c = c+1
    return c

def Consensus(Motifs):
    consensus = ""
    k = len(Motifs[0])
    t = len(Motifs)
    count = CountWithPseudocounts(Motifs)
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus

def CountWithPseudocounts(Motifs):
    count = {}
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(1)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] = count[symbol][j]+1
    return count

def ProfileMostProbablePattern(Text, k, profile):
    n = len(Text)
    m = -1 #This is the bug - need to adjust it to -1 from 0
    x = Text[1:k]
    for i in range(n-k+1):
        Pattern =  Text[i:i+k]
        p = Pr(Pattern, profile)
        if p>m:
            m = p
            x = Pattern
    return x

def Pr(Text, Profile):
    p = 1
    for i in range(len(Text)):
        p1 = Profile[(Text[i])][i]
        p = p*p1
    return p


# GreedyMotifSearch

def Count(Motifs):
    count = {}
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(1)

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
    profile = Count(Motifs)
    for i in 'ACTG':
        for j in range(k):
            profile[i][j] = profile[i][j]/t

    return profile

def Score(Motifs):
    # Insert code here
    score = 0
    k = len(Motifs[0])
    count = Count(Motifs)
    max_symbol = Consensus(Motifs)
    sum1 = 0
    for j in range(k):
        m = 0
        for symbol in "ATCG":
            if count[symbol][j] > m:
                sum1 += count[symbol][j]

    for j in range(k):
        m = 0
        for symbol in "AGTC":
            if count[symbol][j] > m:
                m = count[symbol][j]
        score += m
    return sum1-score

def Pr(Text, Profile):
    p=1
    k = len(Profile["A"])
    for i in range(len(Text)):
        p=p*Profile[Text[i]][i]
    return p

def ProfileMostProbablePattern(text,k,profile):
    p=-1
    result=text[0:k]
    for i in range(len(text)-k+1):
        seq=text[i:i+k]
        pr=Pr(seq,profile)
        if pr>p:
            p=pr
            result=seq

    return result

def GreedyMotifSearch(Dna,k,t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])

    n = len(Dna[0])
    for m in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][m:m+k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))

        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs

    return BestMotifs

Dna = ["GCGCCCCGCCCGGACAGCCATGCGCTAACCCTGGCTTCGATGGCGCCGGCTCAGTTAGGGCCGGAAGTCCCCAATGTGGCAGACCTTTCGCCCCTGGCGGACGAATGACCCCAGTGGCCGGGACTTCAGGCCCTATCGGAGGGCTCCGGCGCGGTGGTCGGATTTGTCTGTGGAGGTTACACCCCAATCGCAAGGATGCATTATGACCAGCGAGCTGAGCCTGGTCGCCACTGGAAAGGGGAGCAACATC", "CCGATCGGCATCACTATCGGTCCTGCGGCCGCCCATAGCGCTATATCCGGCTGGTGAAATCAATTGACAACCTTCGACTTTGAGGTGGCCTACGGCGAGGACAAGCCAGGCAAGCCAGCTGCCTCAACGCGCGCCAGTACGGGTCCATCGACCCGCGGCCCACGGGTCAAACGACCCTAGTGTTCGCTACGACGTGGTCGTACCTTCGGCAGCAGATCAGCAATAGCACCCCGACTCGAGGAGGATCCCG", "ACCGTCGATGTGCCCGGTCGCGCCGCGTCCACCTCGGTCATCGACCCCACGATGAGGACGCCATCGGCCGCGACCAAGCCCCGTGAAACTCTGACGGCGTGCTGGCCGGGCTGCGGCACCTGATCACCTTAGGGCACTTGGGCCACCACAACGGGCCGCCGGTCTCGACAGTGGCCACCACCACACAGGTGACTTCCGGCGGGACGTAAGTCCCTAACGCGTCGTTCCGCACGCGGTTAGCTTTGCTGCC", "GGGTCAGGTATATTTATCGCACACTTGGGCACATGACACACAAGCGCCAGAATCCCGGACCGAACCGAGCACCGTGGGTGGGCAGCCTCCATACAGCGATGACCTGATCGATCATCGGCCAGGGCGCCGGGCTTCCAACCGTGGCCGTCTCAGTACCCAGCCTCATTGACCCTTCGACGCATCCACTGCGCGTAAGTCGGCTCAACCCTTTCAAACCGCTGGATTACCGACCGCAGAAAGGGGGCAGGAC", "GTAGGTCAAACCGGGTGTACATACCCGCTCAATCGCCCAGCACTTCGGGCAGATCACCGGGTTTCCCCGGTATCACCAATACTGCCACCAAACACAGCAGGCGGGAAGGGGCGAAAGTCCCTTATCCGACAATAAAACTTCGCTTGTTCGACGCCCGGTTCACCCGATATGCACGGCGCCCAGCCATTCGTGACCGACGTCCCCAGCCCCAAGGCCGAACGACCCTAGGAGCCACGAGCAATTCACAGCG", "CCGCTGGCGACGCTGTTCGCCGGCAGCGTGCGTGACGACTTCGAGCTGCCCGACTACACCTGGTGACCACCGCCGACGGGCACCTCTCCGCCAGGTAGGCACGGTTTGTCGCCGGCAATGTGACCTTTGGGCGCGGTCTTGAGGACCTTCGGCCCCACCCACGAGGCCGCCGCCGGCCGATCGTATGACGTGCAATGTACGCCATAGGGTGCGTGTTACGGCGATTACCTGAAGGCGGCGGTGGTCCGGA", "GGCCAACTGCACCGCGCTCTTGATGACATCGGTGGTCACCATGGTGTCCGGCATGATCAACCTCCGCTGTTCGATATCACCCCGATCTTTCTGAACGGCGGTTGGCAGACAACAGGGTCAATGGTCCCCAAGTGGATCACCGACGGGCGCGGACAAATGGCCCGCGCTTCGGGGACTTCTGTCCCTAGCCCTGGCCACGATGGGCTGGTCGGATCAAAGGCATCCGTTTCCATCGATTAGGAGGCATCAA", "GTACATGTCCAGAGCGAGCCTCAGCTTCTGCGCAGCGACGGAAACTGCCACACTCAAAGCCTACTGGGCGCACGTGTGGCAACGAGTCGATCCACACGAAATGCCGCCGTTGGGCCGCGGACTAGCCGAATTTTCCGGGTGGTGACACAGCCCACATTTGGCATGGGACTTTCGGCCCTGTCCGCGTCCGTGTCGGCCAGACAAGCTTTGGGCATTGGCCACAATCGGGCCACAATCGAAAGCCGAGCAG", "GGCAGCTGTCGGCAACTGTAAGCCATTTCTGGGACTTTGCTGTGAAAAGCTGGGCGATGGTTGTGGACCTGGACGAGCCACCCGTGCGATAGGTGAGATTCATTCTCGCCCTGACGGGTTGCGTCTGTCATCGGTCGATAAGGACTAACGGCCCTCAGGTGGGGACCAACGCCCCTGGGAGATAGCGGTCCCCGCCAGTAACGTACCGCTGAACCGACGGGATGTATCCGCCCCAGCGAAGGAGACGGCG", "TCAGCACCATGACCGCCTGGCCACCAATCGCCCGTAACAAGCGGGACGTCCGCGACGACGCGTGCGCTAGCGCCGTGGCGGTGACAACGACCAGATATGGTCCGAGCACGCGGGCGAACCTCGTGTTCTGGCCTCGGCCAGTTGTGTAGAGCTCATCGCTGTCATCGAGCGATATCCGACCACTGATCCAAGTCGGGGGCTCTGGGGACCGAAGTCCCCGGGCTCGGAGCTATCGGACCTCACGATCACC"]

# set t equal to the number of strings in Dna and k equal to 15

k = 15
t = 10

# Call GreedyMotifSearch(Dna, k, t) and store the output in a variable called Motifs

Motifs = GreedyMotifSearch(Dna, k, t)

# Print the Motifs variable
print(Motifs)

# Print Score(Motifs)
Print(Score(Motifs))

#Motifs(Profile, Dna) function
def Motifs(Profile, Dna):
    motifs = []
    t = len(Dna)
    k = 4
    for i in range(t):
        motif = ProfileMostProbablePattern(Dna[i], k, Profile)
        motifs.append(motif)
    return motifs
# Insert your ProfileMostProbablePattern(Text, k, Profile) and Pr(Pattern, Profile) functions here.
def Pr(Text, Profile):
    p = 1
    for i in range(len(Text)):
        p = p * Profile[Text[i]][i]
    return p

def ProfileMostProbablePattern(Text, k, Profile):
    p_dict = {}
    for i in range(len(Text)- k +1):
        p = Pr(Text[i: i+k], Profile)
        p_dict[i] = p
    m = max(p_dict.values())
    keys = [k for k,v in p_dict.items() if v == m]
    ind = keys[0]
    return Text[ind: ind +k]


# RandomMotifs (Dna, k, t)
import random

def RandomMotifs(Dna, k, t):

    t = len(Dna)
    l = len(Dna[0])
    RandomMotif =[]
    for i in range(t):
        r = random.randint(0,l-k)
        RandomMotif.append(Dna[i][r:r+k])
    return RandomMotif


# randomized motif search
def randomMotifs(dna,k,t):

    kmm = []
    sc = []
    k = 3
    D = {}
    for i in range(0,len(dna)):
        km = []
        for kk in range(len(dna[i])-k+1):
            km += [dna[i][kk:kk+k]]
        D[i] = km
    for m in range(0,t):
        ran = random.randint(0,len(D[0])-1)
        kmm += [D[m][ran]]

    return kmm

def ProfileWithPseudocounts(Motifs):


    t = len(Motifs)
    k = len(Motifs[0])
    profile = CountWithPseudocounts(Motifs) # output variable
    for symbol in profile:
        for kk in range(0,len(profile[symbol])):
            profile[symbol][kk] = profile[symbol][kk]/(len(Motifs) + 4)

    return profile

def CountWithPseudocounts(Motifs):

    count = {}
    for i in 'ACGT':
        count[i] = []
        for ii in range(len(Motifs[0])):
            count[i].append(1)
    for i in range(len(Motifs)):
        for j in range(len(Motifs[0])):
            symbol = Motifs[i][j]
            count[symbol][j] += 1

    return count

def Score(Motifs):


    count = 0
    L = Consensus(Motifs)
    for i in Motifs:
        for chr1, chr2 in zip(i,L):
            if chr1 != chr2:
                count += 1
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

def Count(Motifs):

    count = {}
    for i in 'ACGT':
        count[i] = []
        for ii in range(len(Motifs[0])):
            count[i].append(0)
    for i in range(len(Motifs)):
        for j in range(len(Motifs[0])):
            symbol = Motifs[i][j]
            count[symbol][j] += 1

    return count

def RandomMotifs(dna,k,t):

    kmm = []
    sc = []
    D = {}
    for i in range(0,len(dna)):
        km = []
        for kk in range(len(dna[i])-k+1):
            km += [dna[i][kk:kk+k]]
        D[i] = km
    for m in range(0,t):
        ran = random.randint(0,len(D[0])-1)
        kmm += [D[m][ran]]

    return kmm

def Motifs(pf,dna):

    k = len(pf['A'])
    D = []
    for i in range(0,len(dna)):
        km = []
        sc = []
        for kk in range(len(dna[i])-k+1):
            km += [dna[i][kk:kk+k]]
        for i in km:
            sc += [Pr(i,pf)]
        D += [km[sc.index(max(sc))]]

    return D

def Pr(Text, Profile):

    p = 1
    for i in range(0,len(Text)):
        p *= Profile[Text[i]][i]

    return p

def RandomizedMotifSearch(Dna, k, t):

    M = RandomMotifs(Dna, k, t)
    BestMotifs = M

    while True:
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile, Dna)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        else:
            return BestMotifs

# Normalize (Probability)
def Normalize(Probabilities):
    # your code here
    sumOfProbabilities = 0
    newlist = {}
    for i in Probabilities:
        sumOfProbabilities += Probabilities.get(i)
    for i in Probabilities:
        newlist[i] = 0
        for j in Probabilities:
            if i == j:
                newlist[i] += Probabilities[j]/sumOfProbabilities
    return newlist

# Weighted Die (Probabilities)

def WeightedDie(Probabilities):
    n = random.uniform(0, 1)
    for p in Probabilities:
        n -= Probabilities[p]
        if n <= 0:
            return p


# ProfileGeneratedString(Text, profile, k)
import random
from operator import itemgetter
# then, copy Pr, Normalize, and WeightedDie below this line
def WeightedDie(d):
    ran = random.uniform(0, 1)
    #print(ran,d)
    tot = 0
    for k, v in sorted(d.items(),key=itemgetter(1)):
        if tot <= ran < v + tot:
            return k
        tot += v
def Normalize(P):
    D = {}
    for k,v in P.items():
        D[k] = P[k]/sum(P.values())
    return D
def Pr(Text, Profile):
    p = 1
    for i in range(0,len(Text)):
        p *= Profile[Text[i]][i]
    return p
# Input:  A string Text, a profile matrix Profile, and an integer k
# Output: ProfileGeneratedString(Text, profile, k)
def ProfileGeneratedString(Text, profile, k):
    # your code here
    n = len(Text)
    probabilities = {}
    for i in range(0,n-k+1):
        probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)

    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)

# GibbsSampler
def GibbsSampler(Dna, k, t, N):
    BestMotifs = [] # output variable
    #randomly select k-mers Motifs = (Motif1, …, Motift) in each string from Dna
    Motifs = RandomMotifs(Dna, k ,t)
    # BestMotifs ← Motifs
    BestMotifs = Motifs.copy()

    #    for j ← 1 to N
    for j in range(N):
        # i ← randomly generated integer between 1 and t
        i = random.randint(0,t-1) # python lists start at 0
        # Profile ← profile matrix formed from all strings in Motifs except for Motifi
        del Motifs[i]
        profile = ProfileWithPseudocounts(Motifs)
        # Motifi ← Profile-randomly generated k-mer in the i-th string
        Motifs.insert(i, ProfileGeneratedString(Dna[i], profile, k))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs
