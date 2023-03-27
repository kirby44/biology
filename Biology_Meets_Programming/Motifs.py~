import random

def Count(Motifs):
    count = {} # initializing the count dictionary
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

def Profile(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}
    for symbol,values in Count(Motifs).items():
        profile[symbol] = [i / t for i in values]
    return profile

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

def Score(Motifs):
    score = 0
    count = Count(Motifs)
    consensus = Consensus(Motifs)
    for i,symbol_c in enumerate(consensus):
        for symbol in "ACGT" :
            if symbol != symbol_c:
                score += count[symbol][i]
    return score

def Pr(Text, Profile):
    p = 1
    for i,symbol in enumerate(Text):
        p *= Profile[symbol][i]
    return p

def ProfileMostProbablePattern(text, k, profile):
    p_max = -1
    probable_kmer = ""
    for i in range(len(text) - k + 1):
        kmer = text[i:i+k]
        p = Pr(kmer, profile)
        if p_max < p:
            p_max = p
            probable_kmer = kmer
    return probable_kmer

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

def CountWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    count = {} # initializing the count dictionary
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

def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}
    for symbol,values in CountWithPseudocounts(Motifs).items():
        profile[symbol] = [i / t for i in values]
    return profile

def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    BestMotifs = [] # output variable
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = ProfileWithPseudocounts(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

def MostProbableMotifs(Profile, Dna, k):
    Motifs = [ProfileMostProbablePattern(Dna[i], k, Profile) for i in range(len(Dna))]
    return Motifs

def RandomMotifs(Dna, k, t):
    import random
    Motifs = []
    for i in range(t):
        r = random.randint(0,len(Dna[i]) - k - 1)
        Motif = Dna[i][r:r+k]
        Motifs.append(Motif)
    return Motifs

def RandomizedMotifSearch(Dna, k, t, n):
    BestMotifs = RandomMotifs(Dna, k, t)
    for i in range(n):
        M = RandomMotifs(Dna, k, t)
        BestMotifs_n = M
        while True:
            Profile = ProfileWithPseudocounts(M)
            M = MostProbableMotifs(Profile, Dna, k)
            if Score(M) < Score(BestMotifs_n):
                BestMotifs_n = M
            else:
                break
        if Score(BestMotifs) > Score(BestMotifs_n):
            BestMotifs = BestMotifs_n
    return BestMotifs

def Normalize(Probabilities):
    s = sum(Probabilities.values())
    for key in Probabilities.keys():
        Probabilities.update({key:Probabilities[key]/s})
    return Probabilities

def WeightedDie(Probabilities):
    r = random.uniform(0,1)
    p = 0
    for kmer,prob in Probabilities.items():
        if p <= r < p + prob:
            return kmer
        else:
            p += prob

def ProfileGeneratedString(Text, profile, k):
    n = len(Text)
    probabilities = {} 
    for i in range(0,n-k+1):
        probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)
    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)

def GibbsSampler(Dna, k, t, N):
    Motifs = RandomMotifs(Dna, k, t)
    BestMotifs = Motifs
    for j in range(N):
        i = random.randint(0,t-1)
        del Motifs[i]
        Profile = ProfileWithPseudocounts(Motifs)
        Motif_i = ProfileGeneratedString(Dna[i], Profile, k)
        Motifs.insert(i, Motif_i)
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

def GibbsSamplerWithRandomStarts(Dna, k, t, N, n):
    Motifs = RandomMotifs(Dna, k, t)
    BestMotifs = Motifs
    for i in range(n):
        M = GibbsSampler(Dna, k, t, N)
        if Score(BestMotifs) > Score(M):
            BestMotifs = M
    return BestMotifs

def MotifEnumeration(Dna, k, d):
    import Replication as Rep
    Patterns = set()
    Neighbor_set_list = []
    t = len(Dna[0])
    for j in range(len(Dna)):
        Neighbors = set()
        for i in range(t - k + 1):
            kmer = Dna[j][i:i+k]
            Neighborhood = Rep.Neighbors(kmer, d)
            Neighbors |= Neighborhood
        Neighbor_set_list.append(Neighbors)
    for Neighbor in Neighbor_set_list[0]:
        for j in range(1,len(Dna)):
            if Neighbor not in Neighbor_set_list[j]:
                break
        else:
            Patterns.add(Neighbor)
    return Patterns

def Entropy(Motifs):
    from math import log2
    entropy = 0
    profile = Profile(Motifs)
    N = len(Motifs[0])
    for i in range(N):
        entropy_i = 0
        for symbol, prs in profile.items():
            pr = prs[i]
            if pr == 0:
                continue
            else:
                entropy_i += -(pr * log2(pr))
        entropy += entropy_i
    return entropy

def DistanceBetweenPatternAndStrings(Pattern, Dna):
    from math import inf
    from Replication import HammingDistance
    k = len(Pattern)
    distance = 0
    for string in Dna:
        HD = inf
        for i in range(len(string) - k + 1):
            kmer = string[i:i+k]
            hd = HammingDistance(Pattern, kmer)
            if HD > hd:
                HD = hd
        distance += HD
    return distance

def AllStrings(k):
    from Replication import Neighbors
    return Neighbors("A" * k, k)

def MedianString(Dna, k):
    from math import inf
    distance = inf
    for kmer in AllStrings(k):
        d = DistanceBetweenPatternAndStrings(kmer, Dna)
        if distance > d:
            distance = d
            Median = kmer
    return Median
            
    
