from Sequence import *
import sys
import os

testdatapath = 'C:/Users/kimur/Documents/testdata/'
ip = 'C:/Users/kimur/Downloads/'
op = 'C:/Users/kimur/Documents/output/'

def test_96_4():
    global ip, op
    ip += sys.argv[1]
    op += sys.argv[1]
    with open(ip, "r") as f:
        lines = f.read().splitlines()
        rna = lines[0]
        O = Translation(rna)
    with open(op, "w") as g:
        g.write(O)
#test_96_4()

def debug_205_5():
    global ip, op
    ip += sys.argv[1]
    inputsdir = ip + 'inputs\\'
    outputsdir = ip + 'outputs\\'
    inputs = os.listdir(inputsdir)
    outputs = os.listdir(outputsdir)
    for input,output in zip(inputs,outputs):
        with open(inputsdir + input, "r") as i, open(outputsdir + output, "r") as o:
            lines = i.read().splitlines()
            kmers = lines[0].split(' ')
            P = ContigGeneration(kmers)
            P.sort()
            output_patterns = o.read().split(' ')
            output_patterns.sort()
            patterns = P
            if patterns == output_patterns:
                print(input," clear")
            else:
                print(patterns)
                print(output_patterns)
                print(input," error")
#debug_205_5()

def week2_q1():
    global ip, op
    ip += sys.argv[1]
    op += sys.argv[1]
    with open(ip, "r") as f:
        lines = f.read().splitlines()
        kmers = []
        for line in lines:
            kmers.append(line)
        P = StringReconstruction(kmers)
        print(P)
    with open(op, "w") as g:
        g.write(P)
#week2_q1()

def week2_q3():
    global ip, op
    ip += sys.argv[1]
    op += sys.argv[1]
    with open(ip, "r") as f:
        lines = f.read().splitlines()
        pairs = []
        k = 3
        d = 1
        for line in lines:
            pairs.append(line)
        print(pairs)
        P = StringReconstructionReadPairs(pairs, k, d)
    with open(op, "w") as g:
        g.write(P)
#week2_q3()

def test_96_7():
    global ip, op
    ip += sys.argv[1]
    op += sys.argv[1]
    with open(ip, "r") as f:
        lines = f.read().splitlines()
        DNA = lines[0]
        peptide = lines[1]
        O = PeptideEncoding(DNA, peptide)
    with open(op, "w") as g:
        for o in O:
            g.write(o + '\n')
#test_96_7()

def test_4912_2():
    global ip, op
    ip += sys.argv[1]
    op += sys.argv[1]
    with open(ip, "r") as f:
        lines = f.read().splitlines()
        peptide = lines[0]
        O = LinearSpectrum(peptide)
        O = [str(i) for i in O]
    with open(op, "w") as g:
        g.write(' '.join(O))
#test_4912_2()

def test_98_4():
    global ip, op
    ip += sys.argv[1]
    op += sys.argv[1]
    with open(ip, "r") as f:
        lines = f.read().splitlines()
        peptide = lines[0]
        O = CyclicSpectrum(peptide)
        O = [str(i) for i in O]
    with open(op, "w") as g:
        g.write(' '.join(O))
#test_98_4()

def test_100_6():
    global ip, op
    ip += sys.argv[1]
    op += sys.argv[1]
    with open(ip, "r") as f:
        lines = f.read().splitlines()
        spectrum = lines[0].split(' ')
        spectrum = [int(i) for i in spectrum]
        O = CyclopeptideSequencing(spectrum)
    with open(op, "w") as g:
        for o in O:
            o = [str(i) for i in o]
            g.write('-'.join(o) + ' ')
#test_100_6()

def test_102_3():
    global ip, op
    ip += sys.argv[1]
    op += sys.argv[1]
    with open(ip, "r") as f:
        lines = f.read().splitlines()
        peptide = lines[0]
        spectrum = lines[1].split(' ')
        spectrum = [int(i) for i in spectrum]
        O = CyclicScore(peptide, spectrum)
    with open(op, "w") as g:
        g.write(str(O))
#test_102_3()

def test_4913_1():
    global ip, op
    ip += sys.argv[1]
    op += sys.argv[1]
    with open(ip, "r") as f:
        lines = f.read().splitlines()
        peptide = lines[0]
        spectrum = lines[1].split(' ')
        spectrum = [int(i) for i in spectrum]
        O = LinearScore(peptide, spectrum)
    with open(op, "w") as g:
        g.write(str(O))
#test_4913_1()

def test_4913_3():
    global ip, op
    ip += sys.argv[1]
    op += sys.argv[1]
    with open(ip, "r") as f:
        lines = f.read().splitlines()
        leaderboard = lines[0].split(' ')
        spectrum = lines[1].split(' ')
        spectrum = [int(i) for i in spectrum]
        N = int(lines[2])
        O = Trim(leaderboard, spectrum, N)
    with open(op, "w") as g:
        g.write(' '.join([str(i) for i in O]))
#test_4913_3()

def test_102_8():
    global ip, op
    ip += sys.argv[1]
    op += sys.argv[1]
    with open(ip, "r") as f:
        lines = f.read().splitlines()
        N = int(lines[0])
        spectrum = lines[1].split(' ')
        spectrum = [int(i) for i in spectrum]
        O = LeaderboardCyclopeptideSequencing(spectrum, N)
        print(O)
    with open(op, "w") as g:
        for o in O:
            g.write('-'.join([str(i) for i in o]) + ' ')
#test_102_8()

def test_104_4():
    global ip, op
    ip += sys.argv[1]
    op += sys.argv[1]
    with open(ip, "r") as f:
        lines = f.read().splitlines()
        spectrum = lines[0].split(' ')
        spectrum = [int(i) for i in spectrum]
        O = Convolution(spectrum)
    with open(op, "w") as g:
        g.write(' '.join([str(i) for i in O]))
#test_104_4()

def test_104_7():
    global ip, op
    ip += sys.argv[1]
    op += sys.argv[1]
    with open(ip, "r") as f:
        lines = f.read().splitlines()
        M = int(lines[0])
        N = int(lines[1])
        spectrum = lines[2].split(' ')
        spectrum = [int(i) for i in spectrum]
        O = ConvolutionCyclopeptideSequencing(spectrum, M, N)
    with open(op, "w") as g:
        for o in O:
            g.write('-'.join([str(i) for i in o]) + ' ')
#test_104_7()

def test_week5():
    global ip, op
    ip += sys.argv[1]
    op += sys.argv[1]
    with open(ip, "r") as f:
        lines = f.read().splitlines()
        length = []
        for line in lines[1:]:
            length.append(int(line.split('\t')[1]))
    sum_length = sum(length)
    long_contigs = len([i for i in length if i >= 1000])
    sum_tops = 0
    for l in length:
        sum_tops += l
        if sum_tops >= sum_length / 2:
            N50 = l
            break
    print('N50 = ',N50)
    print('long contigs = ', long_contigs)
    print('sum length = ', sum_length)
test_week5()



