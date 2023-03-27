from Compare import *
import sys
import os
import numpy as np

testdatapath = 'C:/Users/kimur/Documents/testdata/'
ip = 'C:/Users/kimur/Downloads/'
op = 'C:/Users/kimur/Documents/output/'

def test_261_10():
    global ip, op
    ip += sys.argv[1]
    op += sys.argv[1]
    f = open(ip, "r")
    line = f.readline()[:-1]
    n, m = [int(i) for i in line.split(' ')]
    down = []
    while line:
        line = f.readline()[:-1]
        if line == '-':
            break
        down.append([int(i) for i in line.split(' ')])
    down = np.array(down)
    right = []
    line = f.readline()[:-1]
    while line:
        right.append([int(i) for i in line.split(' ')])
        line = f.readline()[:-1]
    right = np.array(right)
    print(ManhattanTourist(n, m, down, right))
    #with open(op, "w") as g:
    #    g.write(O)
#test_261_10()

def test_245_5():
    global ip, op
    ip += sys.argv[1]
    op += sys.argv[1]
    f = open(ip, "r")
    lines = f.readlines()
    v = lines[0][:-1]
    w = lines[1][:-1]
    backtrack = LCSBackTrack(v, w)
    print(v, w)
    print(backtrack)
    O = OutputLCS(backtrack, v, len(v)-1, len(w)-1)
    print(O)
    print('finish')
#test_245_5()

def test_245_7():
    global ip, op
    ip += sys.argv[1]
    op += sys.argv[1]
    f = open(ip, "r")
    lines = f.readlines()
    source, sink = [int(i) for i in lines[0][:-1].split(' ')]
    edges = []
    for line in lines[1:]:
        edges.append([int(i) for i in line[:-1].split(' ')])
    length, LP = LPBackTrack(edges, source, sink)
    print(length)
    print(LP)
#test_245_7()

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

def test_247_3():
    global ip, op
    ip += sys.argv[1]
    op += sys.argv[1]
    with open(ip, 'r') as f:
        lines = f.read().splitlines()
        match, miss, indel = [int(i) for i in lines[0].split(' ')]
        v = lines[1]
        w = lines[2]
        score, GAv, GAw = GlobalAlignment(v, w, match, miss, indel)
    with open(op, 'w') as o:
        o.write(str(int(score)) +'\n')
        o.write(GAv +'\n')
        o.write(GAw)
#test_247_3()

def debug_247_3():
    global ip, op
    ip += sys.argv[1]
    inputsdir = ip + 'inputs\\'
    outputsdir = ip + 'outputs\\'
    inputs = os.listdir(inputsdir)
    outputs = os.listdir(outputsdir)
    for input,output in zip(inputs,outputs):
        with open(inputsdir + input, "r") as i, open(outputsdir + output, "r") as o:
            lines = i.read().splitlines()
            match, miss, indel = [int(i) for i in lines[0].split(' ')]
            v = lines[1]
            w = lines[2]
            score, GAv, GAw = GlobalAlignment(v, w, match, miss, indel)
            O = o.read().splitlines()
            patterns = str(int(score))
            print(v)
            print(w)
            if patterns == O[0]:
                print(input," clear")
            else:
                print(score, O[0])
                print(input," error")
            print(GAv)
            print(GAw)
            print(O[1])
            print(O[2])
#debug_247_3()

def test_247_10():
    global ip, op
    ip += sys.argv[1]
    op += sys.argv[1]
    with open(ip, 'r') as f:
        lines = f.read().splitlines()
        indel = 5
        v = lines[0]
        w = lines[1]
        score, LAv, LAw = LocalAlignment(v, w, indel)
    with open(op, 'w') as o:
        o.write(str(score) +'\n')
        o.write(LAv +'\n')
        o.write(LAw)
#test_247_10()

def debug_247_10():
    global ip, op
    ip += sys.argv[1]
    inputsdir = ip + 'inputs\\'
    outputsdir = ip + 'outputs\\'
    inputs = os.listdir(inputsdir)
    outputs = os.listdir(outputsdir)
    for input,output in zip(inputs,outputs):
        with open(inputsdir + input, "r") as i, open(outputsdir + output, "r") as o:
            lines = i.read().splitlines()
            indel = 5
            v = lines[0]
            w = lines[1]
            print(v)
            print(w)
            score, LAv, LAw = LocalAlignment(v, w, indel)
            O = o.read().splitlines()
            patterns = str(int(score))
            if patterns == O[0]:
                print(input," clear")
            else:
                print(score, O[0])
                print(input," error")
            print(LAv)
            print(LAw)
            print(O[1])
            print(O[2])
#debug_247_10()

def test_286_4():
    global ip, op
    ip += sys.argv[1]
    op += sys.argv[1]
    with open(ip, 'r') as f:
        lines = f.read().splitlines()
        P = lines[0]
        P = SignToList(P)
        P = GreedySorting(P)
    with open(op, 'w') as o:
        for p in P:
            o.write(p + '\n')
#test_286_4()

def test_289_5():
    global ip, op
    ip += sys.argv[1]
    op += sys.argv[1]
    with open(ip, 'r') as f:
        lines = f.read().splitlines()
        k = int(lines[0])
        P = lines[1]
        Q = lines[2]
        S = Sharedkmers(P, Q, k)
    with open(op, 'w') as o:
        for s in S:
            o.write(str(s) + '\n')
test_289_5()

