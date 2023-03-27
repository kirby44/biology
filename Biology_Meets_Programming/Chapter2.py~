from Sequence import *
import sys
import os

testdatapath = 'C:/Users/kimur/Documents/testdata/'
ip = 'C:/Users/kimur/Downloads/'
op = 'C:/Users/kimur/Documents/output/'

def test_197_3():
    global ip, op
    ip += sys.argv[1]
    op += sys.argv[1]
    with open(ip, "r") as f:
        lines = f.read().splitlines()
        k = int(lines[0])
        Text = lines[1]
        C = Composition(Text, k)
    with open(op, "w") as g:
        g.write('\n'.join(C))
#test_197_3()

def debug_163_4():
    global ip, op
    ip += sys.argv[1]
    #testdir = 'C:\\Users\\kimur\\Documents\\GibbsSampler'
    inputsdir = ip + 'inputs\\'
    outputsdir = ip + 'outputs\\'
    inputs = os.listdir(inputsdir)
    outputs = os.listdir(outputsdir)
    for input,output in zip(inputs,outputs):
        with open(inputsdir + input, "r") as i, open(outputsdir + output, "r") as o:
            lines = i.read().splitlines()
            k = int(lines[0])
            Text = lines[1]
            C = Composition(Text, k)
            output_patterns = o.read().split(" ")
            patterns = C
            if patterns == output_patterns:
                print(input," clear")
            else:
                print(patterns, output_patterns)
                print(input," error")
#debug_163_4()

def debug_198_3():
    global ip, op
    ip += sys.argv[1]
    #testdir = 'C:\\Users\\kimur\\Documents\\GibbsSampler'
    inputsdir = ip + 'inputs\\'
    outputsdir = ip + 'outputs\\'
    inputs = os.listdir(inputsdir)
    outputs = os.listdir(outputsdir)
    for input,output in zip(inputs,outputs):
        with open(inputsdir + input, "r") as i, open(outputsdir + output, "r") as o:
            lines = i.read().splitlines()
            path = lines[0].split(" ")
            G = PathToGenome(path)
            output_patterns = o.read()
            patterns = G
            if patterns == output_patterns:
                print(input," clear")
            else:
                print(patterns, output_patterns)
                print(input," error")
#debug_198_3()

def test_198_3():
    global ip, op
    ip += sys.argv[1]
    op += sys.argv[1]
    with open(ip, "r") as f:
        lines = f.read().splitlines()
        path = lines[0].split(" ")
        G = PathToGenome(path)
    with open(op, "w") as g:
        g.write(G)
#test_198_3()

def test_198_10():
    global ip, op
    ip += sys.argv[1]
    op += sys.argv[1]
    with open(ip, "r") as f:
        lines = f.read().splitlines()
        patterns = lines[0].split(' ')
        O = Overlap(patterns)
    with open(op, "w") as g:
        for key, values in O.items():
            g.write(key + ': ' + ' '.join(values) + '\n')
#test_198_10()

def test_199_6():
    global ip, op
    ip += sys.argv[1]
    op += sys.argv[1]
    with open(ip, "r") as f:
        lines = f.read().splitlines()
        k = int(lines[0])
        Text = lines[1]
        DB = DeBruijnk(Text, k)
    with open(op, "w") as g:
        for key, values in DB.items():
            g.write(key + ': ' + ' '.join(values) + '\n')
#test_199_6()

def test_200_8():
    global ip, op
    ip += sys.argv[1]
    op += sys.argv[1]
    with open(ip, "r") as f:
        lines = f.read().splitlines()
        Patterns = lines[0].split(' ')
        DB = DeBruijnkmers(Patterns)
    with open(op, "w") as g:
        for key, values in DB.items():
            g.write(key + ': ' + ' '.join(values) + '\n')
#test_200_8()

def debug_203_2():
    global ip, op
    ip += sys.argv[1]
    inputsdir = ip + 'inputs\\'
    outputsdir = ip + 'outputs\\'
    inputs = os.listdir(inputsdir)
    outputs = os.listdir(outputsdir)
    for input,output in zip(inputs,outputs):
        with open(inputsdir + input, "r") as i, open(outputsdir + output, "r") as o:
            lines = i.read().splitlines()
            Graph = {}
            for line in lines:
                key, values = line.split(': ')
                key = int(key)
                Graph[key] = [int(j) for j in values.split(' ')]
            output_patterns = [int(j) for j in o.read().split(" ")]
            C = EulerianCycle(Graph)
            patterns = C
            if patterns == output_patterns:
                print(input," clear")
            else:
                print(patterns, output_patterns)
                print(input," error")
#debug_203_2()

def test_203_2():
    global ip, op
    ip += sys.argv[1]
    op += sys.argv[1]
    with open(ip, "r") as f:
        lines = f.read().splitlines()
        Graph = {}
        for line in lines:
            key, values = line.split(': ')
            key = int(key)
            Graph[key] = [int(j) for j in values.split(' ')]
        C = EulerianCycle(Graph)
    with open(op, "w") as g:
        C = [str(j) for j in C]
        g.write(' '.join(C))
#test_203_2()

def test_203_6():
    global ip, op
    ip += sys.argv[1]
    op += sys.argv[1]
    with open(ip, "r") as f:
        lines = f.read().splitlines()
        Graph = {}
        for line in lines:
            key, values = line.split(': ')
            key = int(key)
            Graph[key] = [int(j) for j in values.split(' ')]
        P = EulerianPath(Graph)
    with open(op, "w") as g:
        P = [str(j) for j in P]
        g.write(' '.join(P))
#test_203_6()

def test_203_7():
    global ip, op
    ip += sys.argv[1]
    op += sys.argv[1]
    with open(ip, "r") as f:
        lines = f.read().splitlines()
        Patterns = lines[1].split(' ')
        S = StringReconstruction(Patterns)
    with open(op, "w") as g:
        g.write(S)
#test_203_7()

def test_203_7():
    global ip, op
    ip += sys.argv[1]
    op += sys.argv[1]
    with open(ip, "r") as f:
        lines = f.read().splitlines()
        Patterns = lines[1].split(' ')
        S = StringReconstruction(Patterns)
    with open(op, "w") as g:
        g.write(S)
#test_203_7()

def test_6206_4():
    global ip, op
    ip += sys.argv[1]
    op += sys.argv[1]
    with open(ip, "r") as f:
        lines = f.read().splitlines()
        k, d = [int(i) for i in lines[0].split(' ')]
        GappedPatterns = lines[1].split(' ')
        GappedPatterns = [i.split('|') for i in GappedPatterns]
        S = StringSpelledByGappedPatterns(GappedPatterns, k, d)
    with open(op, "w") as g:
        g.write(S)
#test_6206_4()

def debug_6206_4():
    global ip, op
    ip += sys.argv[1]
    inputsdir = ip + 'inputs\\'
    outputsdir = ip + 'outputs\\'
    inputs = os.listdir(inputsdir)
    outputs = os.listdir(outputsdir)
    for input,output in zip(inputs,outputs):
        with open(inputsdir + input, "r") as i, open(outputsdir + output, "r") as o:
            lines = i.read().splitlines()
            k, d = [int(i) for i in lines[0].split(' ')]
            GappedPatterns = lines[1].split(' ')
            GappedPatterns = [i.split('|') for i in GappedPatterns]
            S = StringSpelledByGappedPatterns(GappedPatterns, k, d)
            output_patterns = o.read()
            patterns = S
            if patterns == output_patterns:
                print(input," clear")
            else:
                print(patterns, output_patterns)
                print(input," error")
#debug_6206_4()

def test_204_16():
    global ip, op
    ip += sys.argv[1]
    op += sys.argv[1]
    with open(ip, "r") as f:
        lines = f.read().splitlines()
        k, d = [int(i) for i in lines[0].split(' ')]
        GappedPatterns = lines[1].split(' ')
        S = StringReconstructionReadPairs(GappedPatterns, k, d)
    with open(op, "w") as g:
        g.write(S)
#test_204_16()

def debug_204_16():
    global ip, op
    ip += sys.argv[1]
    inputsdir = ip + 'inputs\\'
    outputsdir = ip + 'outputs\\'
    inputs = os.listdir(inputsdir)
    outputs = os.listdir(outputsdir)
    for input,output in zip(inputs,outputs):
        with open(inputsdir + input, "r") as i, open(outputsdir + output, "r") as o:
            lines = i.read().splitlines()
            k, d = [int(i) for i in lines[0].split(' ')]
            GappedPatterns = lines[1].split(' ')
            S = StringReconstructionReadPairs(GappedPatterns, k, d)
            output_patterns = o.read()
            patterns = S
            if patterns == output_patterns:
                print(input," clear")
            else:
                print(patterns, output_patterns)
                print(input," error")
#debug_204_16()

def test_205_5():
    global ip, op
    ip += sys.argv[1]
    op += sys.argv[1]
    with open(ip, "r") as f:
        lines = f.read().splitlines()
        kmers = lines[0].split(' ')
        P = ContigGeneration(kmers)
        print(P)
    with open(op, "w") as g:
        g.write(' '.join(P))
#test_205_5()

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
week2_q3()

