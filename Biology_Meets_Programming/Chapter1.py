import os
import Replication as Rep
import Motifs as Mot

def test_2_6():
    Pattern_path = 'C:/Users/kimur/Documents/dataset_2_6.txt'
    with open(Pattern_path,"r") as f:
        lines = f.read().splitlines()
        Text = lines[0]
        Pattern = lines[1]
    print(Rep.PatternCount(Text,Pattern))

def test_2_13():
    testdir = 'C:\\Users\\kimur\\Documents\\FrequentWords\\FrequentWords'
    inputsdir = testdir + '\\inputs\\'
    outputsdir = testdir + '\\outputs\\'
    inputs = os.listdir(inputsdir)
    outputs = os.listdir(outputsdir)
    for input,output in zip(inputs,outputs):
        with open(inputsdir + input, "r") as i, open(outputsdir + output, "r") as o:
            lines = i.read().splitlines()
            Text = lines[0]
            k = int(lines[1])
            output_values = set(o.read().split(" "))
            kmers = set(Rep.BetterFrequentWords(Text, k))
            if kmers == output_values:
                print(input," clear")
            else:
                print(kmers,output_values)
                print(input," error")
#test_2_13()

def test_3_2():
    #testdir = 'C:\\Users\\kimur\\Documents\\ReverseComplement\\ReverseComplement'
    #inputsdir = testdir + '\\inputs\\'
    #outputsdir = testdir + '\\outputs\\'
    #inputs = os.listdir(inputsdir)
    #outputs = os.listdir(outputsdir)
    #for input,output in zip(inputs,outputs):
    #    with open(inputsdir + input, "r") as i, open(outputsdir + output, "r") as o:
    #        lines = i.read().splitlines()
    #        Text = lines[0]
    #        output_values = o.read().split(" ")
    #        reverse = Rep.ReverseComplement(Text)
    #        if reverse == output_values:
    #            print(input," clear")
    #        else:
    #            print(reverse,output_values)
    #            print(input," error")
    Pattern_path = 'C:/Users/kimur/Documents/dataset_3_2.txt'
    with open(Pattern_path,"r") as f:
        lines = f.read().splitlines()
        Text = lines[0]
    print(Rep.ReverseComplement(Text))
#test_3_2()

def test_3_5():
    #Pattern_path = 'C:/Users/kimur/Documents/dataset_3_5.txt'
    #with open(Pattern_path,"r") as f:
    #    lines = f.read().splitlines()
    #    Pattern,Genome = lines
    Pattern_path = 'C:/Users/kimur/Documents/Vibrio_cholerae.txt'
    with open(Pattern_path,"r") as f:
        lines = f.read().splitlines()
        Pattern = 'CTTGATCAT'
        Genome = lines[0]
    print(Rep.PatternMatching(Pattern, Genome))
#test_3_5()

def test_4_5():
    #Pattern_path = 'C:/Users/kimur/Documents/dataset_4_5.txt'
    #with open(Pattern_path,"r") as f:
    #    lines = f.read().splitlines()
    #    Genome = lines[0]
    #    k, L, t = [int(i) for i in lines[1].split()]
    Pattern_path = 'C:/Users/kimur/Documents/E_coli.txt'
    with open(Pattern_path,"r") as f:
        lines = f.read().splitlines()
        Genome = lines[0]
        k, L, t = [9, 500, 3]
    print(Rep.FindClumps(Genome, k, L, t))
#test_4_5()

def week1_quiz():
    #q3
    Text = "TAAACGTGAGAGAAACGTGCTGATTACACTTGTTCGTGTGGTAT"
    print(Rep.FrequentWords(Text, 3))

    #q4
    Pattern = "GATTACA"
    print(Rep.ReverseComplement(Pattern))

    #q5
    Genome = "ATGACTTCGCTGTTACGCGC" 
    Pattern = "CGC"
    print(Rep.PatternMatching(Pattern, Genome))
#week1_quiz()

##Week2
def test_7_8():
    Genome = "GAGCCACCGCGATA"
    print(Rep.SkewArray(Genome))
#test_7_8()

def test_7_10():
    Pattern_path = 'C:/Users/kimur/Documents/dataset_7_10.txt'
    with open(Pattern_path,"r") as f:
        lines = f.read().splitlines()
        Genome = lines[0]
    print(Rep.MinimumSkew(Genome))
#test_7_10()

def test_9_3():
    Pattern_path = 'C:/Users/kimur/Documents/dataset_9_3.txt'
    with open(Pattern_path,"r") as f:
        lines = f.read().splitlines()
        p,q = lines
    print(Rep.HammingDistance(p,q))
#test_9_3()

def test_9_4():
    Pattern_path = 'C:/Users/kimur/Documents/dataset_9_4.txt'
    with open(Pattern_path,"r") as f:
        lines = f.read().splitlines()
        Pattern, Text, d = lines
        d = int(d)
    print(Rep.ApproximatePatternMatching(Pattern, Text, d))
#test_9_4()

def test_9_6():
    Pattern_path = 'C:/Users/kimur/Documents/dataset_9_6.txt'
    with open(Pattern_path,"r") as f:
        lines = f.read().splitlines()
        Pattern, Text, d = lines
        d = int(d)
    #Text = 'AACAAGCTGATAAACATTTAAAGAG'
    #Pattern = 'AAAAA'
    #d = 2
    print(Rep.ApproximatePatternCount(Pattern, Text, d))
#test_9_6()

def test_Neighbors():
    #Pattern = "CAAGTT"
    #d = 2
    Pattern_path = 'C:/Users/kimur/Documents/dataset_3014_4.txt'
    with open(Pattern_path,"r") as f:
        lines = f.read().splitlines()
        Pattern, d = lines
        d = int(d)
    print(Rep.Neighbors(Pattern, d))
#test_Neighbors()

def test_9_9():
    Pattern_path = 'C:/Users/kimur/Documents/dataset_9_9.txt'
    with open(Pattern_path,"r") as f:
        lines = f.read().splitlines()
        Text = lines[0]
        k, d = [int(i) for i in lines[1].split()]
    print(Rep.FrequentWordsWithMismatches(Text, k, d))
#test_9_9()

def test_9_10():
    Pattern_path = 'C:/Users/kimur/Documents/dataset_9_10.txt'
    with open(Pattern_path,"r") as f:
        lines = f.read().splitlines()
        Text = lines[0]
        k, d = [int(i) for i in lines[1].split()]
    print(Rep.FrequentWordsWithMismatchesWithRC(Text, k, d))
#test_9_10()

def week2_quiz():
    #q3
    p = 'CTTGAAGTGGACCTCTAGTTCCTCTACAAAGAACAGGTTGACCTGTCGCGAAG'
    q = 'ATGCCTTACCTAGATGCAATGACGGACGTATTCCTTTTGCCTCAACGGCTCCT'
    print(Rep.HammingDistance(p,q))
    
    #q4
    Genome = 'GCATACACTTCCCAGTAGGTACTG'
    array = Rep.SkewArray(Genome)
    max_Skew = max(array)
    positions = [i for i in range(len(array)) if array[i] == max_Skew]
    print(array)
    print(positions)

    #q5
    Text = 'CGTGACAGTGTATGGGCATCTTT'
    Pattern = 'TGT'
    d = 1
    print(Rep.ApproximatePatternCount(Pattern, Text, d))

    #q6
    Pattern = 'CCAGTCAATG'
    d = 1
    print(len(Rep.Neighbors(Pattern, d)))
#week2_quiz()

#week3
def debug_156_8():
    testdir = 'C:\\Users\\kimur\\Documents\\MotifEnumeration\\MotifEnumeration'
    inputsdir = testdir + '\\inputs\\'
    outputsdir = testdir + '\\outputs\\'
    inputs = os.listdir(inputsdir)
    outputs = os.listdir(outputsdir)
    for input,output in zip(inputs,outputs):
        with open(inputsdir + input, "r") as i, open(outputsdir + output, "r") as o:
            lines = i.read().splitlines()
            k, d = [int(i) for i in lines[0].split()]
            Dna = lines[1].split()
            output_patterns = set(o.read().split(" "))
            patterns = Mot.MotifEnumeration(Dna, k, d)
            if patterns == output_patterns:
                print(input," clear")
            else:
                print(patterns, output_patterns)
                print(input," error")
#debug_156_8()

def test_156_8():
    Pattern_path = 'C:/Users/kimur/Documents/dataset_156_8.txt'
    with open(Pattern_path,"r") as f:
        lines = f.read().splitlines()
        k, d = [int(i) for i in lines[0].split()]
        Dna = lines[1].split()
    print(Mot.MotifEnumeration(Dna, k, d))
#test_156_8()

def test_157_8():
    Motifs = [
    "TCGGGGGTTTTT",
    "CCGGTGACTTAC",
    "ACGGGGATTTTC",
    "TTGGGGACTTTT",
    "AAGGGGACTTCC",
    "TTGGGGACTTCC",
    "TCGGGGATTCAT",
    "TCGGGGATTCCT",
    "TAGGGGAACTAC",
    "TCGGGTATAACC"
    ]
    print(Mot.Entropy(Motifs))
#test_157_8()

def test_5164():
    Pattern_path = 'C:/Users/kimur/Documents/dataset_5164_1.txt'
    with open(Pattern_path,"r") as f:
        lines = f.read().splitlines()
        Pattern = lines[0]
        Dna = lines[1].split()
    print(Mot.DistanceBetweenPatternAndStrings(Pattern, Dna))
#test_5164()

def test_158_9():
    Pattern_path = 'C:/Users/kimur/Documents/dataset_158_9.txt'
    with open(Pattern_path,"r") as f:
        lines = f.read().splitlines()
        k = int(lines[0])
        Dna = lines[1].split()
    print(Mot.MedianString(Dna, k))
#test_158_9()

def test_159_3():
    Pattern_path = 'C:/Users/kimur/Documents/dataset_159_3.txt'
    with open(Pattern_path,"r") as f:
        lines = f.read().splitlines()
        Text = lines[0]
        k = int(lines[1])
        P = lines[2:6]
        profile = {}
        for i, symbol in enumerate(["A","C","G","T"]):
            profile[symbol] = [float(j) for j in P[i].split(' ')]
    print(Mot.ProfileMostProbablePattern(Text, k, profile))
#test_159_3()

def test_159_5():
    Pattern_path = 'C:/Users/kimur/Documents/dataset_159_5.txt'
    with open(Pattern_path,"r") as f:
        lines = f.read().splitlines()
        k, t = [int(i) for i in lines[0].split(' ')]
        Dna = lines[1].split(' ')
    print(Mot.GreedyMotifSearch(Dna, k, t))
#test_159_5()

def test_160_9():
    Pattern_path = 'C:/Users/kimur/Documents/dataset_160_9.txt'
    with open(Pattern_path,"r") as f:
        lines = f.read().splitlines()
        k, t = [int(i) for i in lines[0].split(' ')]
        Dna = lines[1].split(' ')
    print(Mot.GreedyMotifSearchWithPseudocounts(Dna, k, t))
#test_160_9()

#quiz3
def week3_quiz():
    Dna = [
            'CTCGATGAGTAGGAAAGTAGTTTCACTGGGCGAACCACCCCGGCGCTAATCCTAGTGCCC',
            'GCAATCCTACCCGAGGCCACATATCAGTAGGAACTAGAACCACCACGGGTGGCTAGTTTC',
            'GGTGTTGAACCACGGGGTTAGTTTCATCTATTGTAGGAATCGGCTTCAAATCCTACACAG'
            ]
    k = 7
    print(Mot.MedianString(Dna, k))
#week3_quiz()

def test_161_5():
    Pattern_path = 'C:/Users/kimur/Documents/dataset_161_5.txt'
    with open(Pattern_path,"r") as f:
        lines = f.read().splitlines()
        k, t = [int(i) for i in lines[0].split(' ')]
        Dna = lines[1].split(' ')
        n = 1000
        M = Mot.RandomizedMotifSearch(Dna, k, t, n)
        print(M)
#test_161_5()

def debug_161_5():
    testdir = 'C:\\Users\\kimur\\Documents\\RandomizedMotifSearch'
    inputsdir = testdir + '\\inputs\\'
    outputsdir = testdir + '\\outputs\\'
    inputs = os.listdir(inputsdir)
    outputs = os.listdir(outputsdir)
    for input,output in zip(inputs,outputs):
        with open(inputsdir + input, "r") as i, open(outputsdir + output, "r") as o:
            lines = i.read().splitlines()
            k, t = [int(i) for i in lines[0].split(' ')]
            Dna = lines[1].split(' ')
            n = 1000
            output_patterns = set(o.read().split(" "))
            patterns = set(Mot.RandomizedMotifSearch(Dna, k, t, n))
            if patterns == output_patterns:
                print(input," clear")
            else:
                print(patterns, output_patterns)
                print(input," error")
#debug_161_5()

def test_163_4():
    Pattern_path = 'C:/Users/kimur/Documents/dataset_163_4.txt'
    with open(Pattern_path,"r") as f:
        lines = f.read().splitlines()
        k, t, N = [int(i) for i in lines[0].split(' ')]
        Dna = lines[1].split(' ')
        n = 20
        M = Mot.GibbsSamplerWithRandomStarts(Dna, k, t, N, n)
        print(M)
test_163_4()

def debug_163_4():
    testdir = 'C:\\Users\\kimur\\Documents\\GibbsSampler'
    inputsdir = testdir + '\\inputs\\'
    outputsdir = testdir + '\\outputs\\'
    inputs = os.listdir(inputsdir)
    outputs = os.listdir(outputsdir)
    for input,output in zip(inputs,outputs):
        with open(inputsdir + input, "r") as i, open(outputsdir + output, "r") as o:
            lines = i.read().splitlines()
            k, t, N = [int(i) for i in lines[0].split(' ')]
            n = 20
            Dna = lines[1].split(' ')
            output_patterns = set(o.read().split(" "))
            M = Mot.GibbsSamplerWithRandomStarts(Dna, k, t, N, n)
            patterns = set(M)
            if patterns == output_patterns:
                print(input," clear")
            else:
                print(patterns, output_patterns)
                print(input," error")
#debug_163_4()

