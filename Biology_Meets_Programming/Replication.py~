def PatternCount(Text,Pattern):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count

def FrequencyMap(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern] = 0
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern] += 1
    return freq

def FindClumps(Text, k, L, t):
    Patterns = set()
    n = len(Text)
    for i in range(n - L):
        Window = Text[i:i+L]
        freqMap = FrequencyMap(Window, k)
        for kmer,freq in freqMap.items():
            if freq >= t:
                Patterns.add(kmer)
    return Patterns
        
def BetterFrequentWords(Text, k):
    frequentPatterns = []
    freqMap = FrequencyMap(Text, k)
    max_freq = max(freqMap.values())
    for kmer,freq in freqMap.items():
        if freq == max_freq:
            frequentPatterns.append(kmer)
    return frequentPatterns

def FrequentWords(Text, k):
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        if freq[key] == m:
            words.append(key)
    return words

def Reverse(Pattern):
    rev = ''
    for char in Pattern:
        rev = char + rev
    return rev

def Complement(Pattern):
    complement_char = {"A":"T","T":"A","G":"C","C":"G"}
    com = ''
    for char in Pattern:
        com += complement_char[char]
    return com

def Transcription(DNA):
    RNA = DNA.replace('T', 'U')
    return RNA
        
def ReverseComplement(Pattern):   
    Pattern = Reverse(Pattern)
    Pattern = Complement(Pattern)
    return Pattern

def PatternMatching(Pattern, Genome):
    positions = [] # output variable
    k = len(Pattern)
    for i in range(len(Genome) - k +1):
        if Genome[i:i + k] == Pattern:
            positions.append(i)
    return positions

def SymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        array[i] = PatternCount(symbol, ExtendedGenome[i:i+(n//2)])
    return array

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
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array

def SkewArray(Genome):
    array = [0]
    for i in range(len(Genome)):
        chr = Genome[i]
        if chr in ["A" ,"T"]:
            array.append(array[i])
        if chr == "G":
            num = array[i]+1
            array.append(num)
        if chr == "C":
            num = array[i]-1
            array.append(num)
    return array

def MinimumSkew(Genome):
    array = SkewArray(Genome)
    min_Skew = min(array)
    positions = [i for i in range(len(array)) if array[i] == min_Skew]
    return positions

def HammingDistance(p, q):
    count = 0
    for i,j in zip(p,q):
        if i != j:
            count += 1
    return count

def ApproximatePatternMatching(Pattern, Text, d):
    positions = []
    for i in range(len(Text) - len(Pattern) + 1):
        if HammingDistance(Pattern,Text[i:i+len(Pattern)]) <= d:
            positions.append(i)
    return positions

def ApproximatePatternCount(Pattern, Text, d):
    positions = []
    for i in range(len(Text) - len(Pattern) + 1):
        if HammingDistance(Pattern,Text[i:i+len(Pattern)]) <= d:
            positions.append(i)
    return len(positions)

def Neighbors(Pattern, d):
    if d == 0:
        return {Pattern}
    if len(Pattern) == 1:
        return {"A", "C", "G", "T"}
    Neighborhood = set()
    Suffix = Pattern[1:]
    SuffixNeighbors = Neighbors(Suffix, d)
    for SuffixNeighbor in SuffixNeighbors:
        if HammingDistance(Suffix, SuffixNeighbor) < d:
            for x in ["A","G","T","C"]:
                Neighborhood.add(x + SuffixNeighbor)
        else:
                Neighborhood.add(Pattern[0] + SuffixNeighbor)
    return Neighborhood

def FrequentWordsWithMismatches(Text, k, d):
    Patterns = []
    freqMap = {}
    n = len(Text)
    for i in range(n-k):
        Pattern = Text[i:i+k]
        neighborhood = Neighbors(Pattern, d)
        for neighbor in neighborhood:
            if not neighbor in freqMap:
                freqMap[neighbor] = 1
            else:
                freqMap[neighbor] += 1
    m = max(freqMap.values())
    for Pattern, freq in freqMap.items():
        if freq == m:
            Patterns.append(Pattern)
    return Patterns

def FrequentWordsWithMismatchesWithRC(Text, k, d):
    Patterns = []
    freqMap = {}
    freqMapwithRC = {}
    n = len(Text)
    for i in range(n-k):
        Pattern = Text[i:i+k]
        neighborhood = Neighbors(Pattern, d)
        for neighbor in neighborhood:
            if not neighbor in freqMap:
                freqMap[neighbor] = 1
            else:
                freqMap[neighbor] += 1
    for Pattern, freq in freqMap.items():
        if ReverseComplement(Pattern) == Pattern:
            freqMapwithRC[Pattern] = freqMap[Pattern]
            continue
        elif not ReverseComplement(Pattern) in freqMap:
            freqMapwithRC[Pattern] = freqMap[Pattern]
            continue
        else:
            freqMapwithRC[Pattern] = freqMap[Pattern] + freqMap[ReverseComplement(Pattern)]
    m = max(freqMapwithRC.values())
    for Pattern, freq in freqMapwithRC.items():
        if freq == m:
            Patterns.append(Pattern)
    return Patterns
