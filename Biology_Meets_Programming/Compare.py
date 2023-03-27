import sys
import math
import numpy as np
import pandas as pd
import itertools

def MinNumCoins(m):
    MinNumCoins = [0, 1, 2, 3, 1]
    denominations = [1, 4, 5]
    for mi in range(5, m + 1):
        MinNumCoins.append(min([MinNumCoins[mi - d] + 1 for d in denominations]))
    return MinNumCoins

def DPChange(money, Coins):
    MinNumCoins = [0]
    for m in range(1, money + 1):
        MinNumCoins.append(math.inf)
        for C in Coins:
            if m >= C:
                NumC = MinNumCoins[m - C] + 1
                if MinNumCoins[m] > NumC:
                    MinNumCoins[m] = NumC
    return MinNumCoins[money]

def ManhattanTourist(n, m, Down, Right):
    S = np.zeros((n+1, m+1))
    for i in range(n):
        S[i+1, 0] = S[i, 0] + Down[i, 0]
    for j in range(m):
        S[0, j+1] = S[0, j] + Right[0, j]
    for i in range(n):
        for j in range(m):
            D = S[i, j+1] + Down[i, j+1]
            R = S[i+1, j] + Right[i+1, j]
            S[i+1, j+1] = max(D, R)
    return S[n, m]

def LCSBackTrack(v, w):
    Backtrack = np.empty((len(v), len(w)), dtype = 'object')
    S = np.zeros((len(v)+1, len(w)+1))
    for i in range(len(v)):
        for j in range(len(w)):
            match = 0
            if v[i] == w[j]:
                match = 1
            S[i+1, j+1] = max(S[i, j+1], S[i+1, j], S[i, j] + match)
            if S[i+1, j+1] == S[i, j+1]:
                Backtrack[i, j] = 'd'
            elif S[i+1, j+1] == S[i+1, j]:
                Backtrack[i, j] = 'r'
            elif S[i+1, j+1] == S[i, j] + match:
                Backtrack[i, j] = 'x'
    return Backtrack

def OutputLCS(backtrack, v, i, j):
    LCS = ''
    while i > 0 and j > 0:
        B = backtrack[i, j]
        if B == 'd':
            i -= 1
        elif B == 'r':
            j -= 1
        elif B == 'x':
            LCS = v[i] + LCS
            i -= 1
            j -= 1
    else:
        if B == 'x':
            LCS = v[i] + LCS
    return LCS

    #sys.setrecursionlimit(100000000)
    #if B == 'd':
    #    return OutputLCS(backtrack, v, i-1, j)
    #elif B == 'r':
    #    return OutputLCS(backtrack, v, i, j-1)
    #elif B == 'x':
    #    print(v[i])
    #    return OutputLCS(backtrack, v, i-1, j-1) + v[i]

def LPBackTrack(edges, source, sink):
    print('source -sink = ', source, '-', sink)
    print('nedges', len(edges))
    innodes = [edge[0] for edge in edges]
    outnodes = [edge[1] for edge in edges]
    nodes = set(innodes + outnodes)
    Backtrack = dict.fromkeys(nodes)
    S = dict.fromkeys(nodes)
    S[source] = 0
    while S[sink] == None:
        for node in nodes:
            if node == source:
                continue
            if not S[node] == None:
                continue
            predecessors = []
            S_ps = []
            is_ready = True
            for edge in edges:
                if edge[1] == node:
                    if S[edge[0]] == None:
                        is_ready = False
                        break
                    S_ps.append(S[edge[0]] + edge[2])
                    predecessors.append(edge[0])
            if is_ready:
                S[node] = max(S_ps)
            for edge in edges:
                if edge[1] == node:
                    if S[edge[0]] + edge[2] == S[node]:
                        Backtrack[node] = edge[0]
    LP = str(sink)
    node = sink
    while Backtrack[node] != source:
        LP = str(Backtrack[node]) + ' ' + LP
        node = Backtrack[node]
    LP = str(source) + ' ' + LP
    return S[sink], LP

def GlobalAlignment(v, w, match_score, mismatch_score, indel_score):
    backtrack = np.empty((len(v)+1, len(w)+1), dtype = 'object')
    backtrack[1:,0] = np.array(['d'] * len(v))
    backtrack[0,1:] = np.array(['r'] * len(w))
    S = np.zeros((len(v)+1, len(w)+1))
    S[:,0] = np.array([i * (-indel_score) for i in range(len(v)+1)])
    S[0,:] = np.array([i * (-indel_score) for i in range(len(w)+1)])
    for i in range(len(v)):
        for j in range(len(w)):
            match = 0
            if v[i] == w[j]:
                match = match_score
            else:
                match = - mismatch_score
            Si = S[i, j+1] - indel_score
            Sj = S[i+1, j] - indel_score
            Sij = S[i, j] + match
            S[i+1, j+1] = max(Si, Sj, Sij)
            if S[i+1, j+1] == Si:
                if i > 0:
                    backtrack[i+1, j+1] = 'd'
            elif S[i+1, j+1] == Sj:
                backtrack[i+1, j+1] = 'r'
            elif S[i+1, j+1] == Sij:
                    backtrack[i+1, j+1] = 'x'
    print(backtrack)
    print(S)
    GAv = ''
    GAw = ''
    i = len(v)-1
    j = len(w)-1
    while i >= 0 or j >= 0:
        B = backtrack[i+1, j+1]
        if B == 'd':
            GAv = v[i] + GAv
            GAw = '-' + GAw
            i -= 1
        elif B == 'r':
            GAv = '-' + GAv
            GAw = w[j] + GAw
            j -= 1
        elif B == 'x':
            GAv = v[i] + GAv
            GAw = w[j] + GAw
            i -= 1
            j -= 1
    return int(S[len(v),len(w)]), GAv, GAw

def LocalAlignment(v, w, indel_score):
    PAM_path = 'C:/Users/kimur/Documents/src/Biology_Meets_Programming/PAM250.txt'
    PAM = pd.read_csv(PAM_path, delim_whitespace=True)
    backtrack = np.empty((len(v)+1, len(w)+1), dtype = 'object')
    backtrack[0,0] = 's'
    backtrack[1:,0] = np.array(['d'] * len(v))
    backtrack[0,1:] = np.array(['r'] * len(w))
    S = np.zeros((len(v)+1, len(w)+1))
    for i in range(len(v)):
        for j in range(len(w)):
            match = PAM[v[i]][w[j]]
            Si = S[i, j+1] - indel_score
            Sj = S[i+1, j] - indel_score
            Sij = S[i, j] + match
            S[i+1, j+1] = max(Si, Sj, Sij, 0)
            if S[i+1, j+1] == Si:
                if i > 0:
                    backtrack[i+1, j+1] = 'd'
            elif S[i+1, j+1] == Sj:
                backtrack[i+1, j+1] = 'r'
            elif S[i+1, j+1] == Sij:
                backtrack[i+1, j+1] = 'x'
            elif int(S[i+1, j+1]) == 0:
                backtrack[i+1, j+1] = 's'
    print(backtrack)
    print(S)
    LAv = ''
    LAw = ''
    i,j = [a[0]-1 for a in np.where(S == S.max())] 
    while backtrack[i+1, j+1] != 's':
        B = backtrack[i+1, j+1]
        if B == 'd':
            LAv = v[i] + LAv
            LAw = '-' + LAw
            i -= 1
        elif B == 'r':
            LAv = '-' + LAv
            LAw = w[j] + LAw
            j -= 1
        elif B == 'x':
            LAv = v[i] + LAv
            LAw = w[j] + LAw
            i -= 1
            j -= 1
    return int(S[len(v),len(w)]), LAv, LAw

def editDistance(x, y):
    backtrack = np.empty((len(x)+1, len(y)+1), dtype = 'object')
    backtrack[0,0] = 's'
    backtrack[1:,0] = np.array(['v'] * len(x))
    backtrack[0,1:] = np.array(['h'] * len(y))
    # Create distance matrix
    D = []
    for i in range(len(x)+1):
        D.append([0]*(len(y)+1))
    # Initialize first row and column of matrix
    for i in range(len(x)+1):
        D[i][0] = i
    for i in range(len(y)+1):
        D[0][i] = i
    # Fill in the rest of the matrix
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            distHor = D[i][j-1] + 1
            distVer = D[i-1][j] + 1
            if x[i-1] == y[j-1]:
                distDiag = D[i-1][j-1]
            else:
                distDiag = D[i-1][j-1] + 1
            D[i][j] = min(distHor, distVer, distDiag)
            if D[i][j] == distHor:
                backtrack[i, j] = 'h'
            elif D[i][j] == distVer:
                backtrack[i, j] = 'v'
            elif D[i][j] == distDiag:
                backtrack[i, j] = 'd'
    Fx = ''
    Fy = ''
    i = len(x)-1
    j = len(y)-1
    print(backtrack)
    print(D)
    print(D[-1][-1])
    B = ''
    while B != 's':
        B = backtrack[i+1, j+1]
        print(i,j,B)
        if B == 'h':
            Fx = x[i] + Fx
            Fy = '-' + Fy
            j -= 1
        elif B == 'v':
            Fx = '-' + Fx
            Fy = y[j] + Fy
            i -= 1
        elif B == 'd':
            Fx = x[i] + Fx
            Fy = y[j] + Fy
            i -= 1
            j -= 1
    # Edit distance is the value in the bottom right corner of the matrix
    return D[-1][-1], Fx, Fy

def FittingAlignment(v, w, indel_score):
    BLOSUM_path = 'C:/Users/kimur/Documents/src/Biology_Meets_Programming/BLOSUM62.txt'
    BLOSUM = pd.read_csv(BLOSUM_path, delim_whitespace=True)
    backtrack = np.empty((len(v)+1, len(w)+1), dtype = 'object')
    backtrack[0,0] = 's'
    backtrack[1:,0] = np.array(['s'] * len(v))
    backtrack[0,1:] = np.array(['s'] * len(w))
    S = np.zeros((len(v)+1, len(w)+1))
    for i in range(len(v)):
        for j in range(len(w)):
            match = BLOSUM[v[i]][w[j]]
            Si = S[i, j+1] - indel_score
            Sj = S[i+1, j] - indel_score
            Sij = S[i, j] + match
            S[i+1, j+1] = max(Si, Sj, Sij, 0)
            if S[i+1, j+1] == Si:
                if i > 0:
                    backtrack[i+1, j+1] = 'd'
            elif S[i+1, j+1] == Sj:
                backtrack[i+1, j+1] = 'r'
            elif S[i+1, j+1] == Sij:
                backtrack[i+1, j+1] = 'x'
            elif int(S[i+1, j+1]) == 0:
                backtrack[i+1, j+1] = 's'
    print(backtrack)
    print(S)
    LAv = ''
    LAw = ''
    lastcolumn = S[:,len(w)]
    score = lastcolumn.max()
    i = np.where(lastcolumn == score)[0][0] - 1
    print(i)
    j = len(w) - 1
    while backtrack[i+1, j+1] != 's':
        B = backtrack[i+1, j+1]
        if B == 'd':
            LAv = v[i] + LAv
            LAw = '-' + LAw
            i -= 1
        elif B == 'r':
            LAv = '-' + LAv
            LAw = w[j] + LAw
            j -= 1
        elif B == 'x':
            LAv = v[i] + LAv
            LAw = w[j] + LAw
            i -= 1
            j -= 1
    return int(score), LAv, LAw


def OverlapAlignment(v, w, match_score, mismatch_score, indel_score):
    backtrack = np.empty((len(v)+1, len(w)+1), dtype = 'object')
    backtrack[0,0] = 's'
    backtrack[1:,0] = np.array(['s'] * len(v))
    backtrack[0,1:] = np.array(['r'] * len(w))
    S = np.zeros((len(v)+1, len(w)+1))
    S[0,:] = np.array([i * -indel_score for i in range(len(w)+1)])
    for i in range(len(v)):
        for j in range(len(w)):
            if v[i] == w[j]:
                match = match_score
            else:
                match = -mismatch_score
            Si = S[i, j+1] - indel_score
            Sj = S[i+1, j] - indel_score
            Sij = S[i, j] + match
            S[i+1, j+1] = max(Si, Sj, Sij)
            if S[i+1, j+1] == Si:
                if i > 0:
                    backtrack[i+1, j+1] = 'd'
            elif S[i+1, j+1] == Sj:
                backtrack[i+1, j+1] = 'r'
            elif S[i+1, j+1] == Sij:
                backtrack[i+1, j+1] = 'x'
            #elif j == 0 and int(S[i+1, j+1]) == 0:
            #    backtrack[i+1, j+1] = 's'
    print(backtrack)
    print(S)
    LAv = ''
    LAw = ''
    lastrow = S[len(v),:]
    score = lastrow.max()
    i = len(v) - 1
    j = np.where(lastrow == score)[0][-1] - 1
    while backtrack[i+1, j+1] != 's':
        B = backtrack[i+1, j+1]
        if B == 'd':
            LAv = v[i] + LAv
            LAw = '-' + LAw
            i -= 1
        elif B == 'r':
            LAv = '-' + LAv
            LAw = w[j] + LAw
            j -= 1
        elif B == 'x':
            LAv = v[i] + LAv
            LAw = w[j] + LAw
            i -= 1
            j -= 1
    return int(score), LAv, LAw

def AffineAlignment(v, w, match_score, mismatch_score, gap_opening_score, gap_extension_score):
    backtrack = np.empty((len(v)+1, len(w)+1), dtype = 'object')
    backtrack[0,0] = 's'
    backtrack[1:,0] = np.array(['d'] * len(v))
    backtrack[0,1:] = np.array(['r'] * len(w))
    Sl = np.zeros((len(v)+1, len(w)+1))
    Sl[0,:] = np.full(len(w)+1, -np.inf)
    Sl[:,0] = np.full(len(v)+1, -np.inf)
    Sm = np.zeros((len(v)+1, len(w)+1))
    Sm[0,1:] = np.array([-gap_extension_score * i - gap_opening_score for i in range(len(w))])
    Sm[1:,0] = np.array([-gap_extension_score * i - gap_opening_score for i in range(len(v))])
    Su = np.zeros((len(v)+1, len(w)+1))
    Su[0,:] = np.full(len(w)+1, -np.inf)
    Su[:,0] = np.full(len(v)+1, -np.inf)
    for i in range(len(v)):
        for j in range(len(w)):
            if v[i] == w[j]:
                match = match_score
            else:
                match = -mismatch_score
            Sl_open = Sm[i, j+1] - gap_opening_score 
            Sl_exte = Sl[i, j+1] - gap_extension_score
            Sl[i+1, j+1] = max(Sl_open, Sl_exte)
            Su_open = Sm[i+1, j] - gap_opening_score
            Su_exte = Su[i+1, j] - gap_extension_score
            Su[i+1, j+1] = max(Su_open, Su_exte)
            Sm_lowe = Sl[i+1, j+1]
            Sm_matc = Sm[i, j] + match
            Sm_uppe = Su[i+1, j+1]
            Sm[i+1, j+1] = s = max(Sm_lowe, Sm_matc, Sm_uppe)
            if s == Sm_uppe:
                backtrack[i+1, j+1] = 'r'
            elif s == Sm_lowe:
                backtrack[i+1, j+1] = 'd'
            elif s == Sm_matc:
                backtrack[i+1, j+1] = 'x'
    print(backtrack)
    print(Sl)
    print(Sm)
    print(Su)
    LAv = ''
    LAw = ''
    score = Sm[-1,-1]
    i = len(v) - 1
    j = len(w) - 1
    while backtrack[i+1, j+1] != 's':
        B = backtrack[i+1, j+1]
        if B == 'd':
            LAv = v[i] + LAv
            LAw = '-' + LAw
            i -= 1
        elif B == 'r':
            LAv = '-' + LAv
            LAw = w[j] + LAw
            j -= 1
        elif B == 'x':
            LAv = v[i] + LAv
            LAw = w[j] + LAw
            i -= 1
            j -= 1
    return int(score), LAv, LAw

def Si_GlobalAlignment(v, w, match_score, mismatch_score, indel_score):
    S = np.zeros((len(v)+1, len(w)+1))
    S[:,0] = np.array([i * (-indel_score) for i in range(len(v)+1)])
    S[0,:] = np.array([i * (-indel_score) for i in range(len(w)+1)])
    for i in range(len(v)):
        for j in range(len(w)):
            match = 0
            if v[i] == w[j]:
                match = match_score
            else:
                match = - mismatch_score
            Si = S[i, j+1] - indel_score
            Sj = S[i+1, j] - indel_score
            Sij = S[i, j] + match
            S[i+1, j+1] = max(Si, Sj, Sij)
    return S[:,-1]

def MiddleEdge(v, w, match_score, mismatch_score, indel_score):
    j_middle =  len(w) // 2
    w_left = w[:j_middle]
    w_right = w[j_middle+1:]
    Si_FromSource = Si_GlobalAlignment(v, w_left, match_score, mismatch_score, indel_score)
    Sn_minus_i_ToSink = Si_GlobalAlignment(v[::-1], w_right[::-1], match_score, mismatch_score, indel_score)
    Si_ToSink = np.flip(Sn_minus_i_ToSink)
    Length = -math.inf
    for i in range(len(v)):
        print(i,j_middle)
        if v[i] == w[j_middle]:
            match = match_score
        else:
            match = -mismatch_score
        Length_i_right = Si_FromSource[i] -indel_score + Si_ToSink[i] 
        if Length < Length_i_right:
            Length = Length_i_right
            i_middle_from = i
            i_middle_to = i
            edge = 'r'
        if i == len(v):
            continue
        Length_i_diag = Si_FromSource[i] + match + Si_ToSink[i+1] 
        if Length < Length_i_diag:
            Length = Length_i_diag
            i_middle_from = i
            i_middle_to = i+1
            edge = 'x'
    node_from = (i_middle_from, j_middle)
    node_to = (i_middle_to, j_middle+1)
    print(node_from, node_to)

def MultipleAlignment(v, w, u):
    S = np.zeros((len(v)+1, len(w)+1, len(u)+1))
    backtrack = np.full((len(v)+1, len(w)+1, len(u)+1), '000') 
    backtrack[0, 0, 0] = 's'
    for i in range(len(v)):
        for j in range(len(w)):
            for k in range(len(u)):
                I = '0'
                J = '0'
                K = '0'
                if i == 0:
                    I = '1'
                if j == 0:
                    J = '1'
                if k == 0:
                    K = '1'
                backtrack[i, j, k] = I + J + K
    for i in range(len(v)):
        for j in range(len(w)):
            for k in range(len(u)):
                candidates = []
                for I, J, K in itertools.product((0, 1), repeat=3):
                    if (I, J, K) != (0, 0, 0) and (I, J, K) != (1, 1, 1): 
                        candidates.append(S[i+I, j+J, k+K])
                if v[i] == w[j] == u[k]:
                    candidates.append(S[i, j, k] + 1)
                else:
                    candidates.append(S[i, j, k])
                S[i+1, j+1, k+1] = s = max(candidates)
                if S[i, j, k] + 1 == s:
                    backtrack[i+1, j+1, k+1] = '000'
                for I, J, K in itertools.product((0, 1), repeat=3):
                    if (I, J, K) != (1, 1, 1): 
                        if S[i+I, j+J, k+K] == s:
                            backtrack[i+1, j+1, k+1] = str(I) + str(J) + str(K)
                            break
    MAv = ''
    MAw = ''
    MAu = ''
    i = len(v)-1
    j = len(w)-1
    k = len(u)-1
    print(S)
    print(backtrack)
    B = backtrack[i+1, j+1, k+1]
    while B != '111':
        if B[0] == '0':
            MAv = v[i] + MAv
            i -= 1
        elif B[0] == '1':
            MAv = '-' + MAv
        if B[1] == '0':
            MAw = w[j] + MAw
            j -= 1
        elif B[1] == '1':
            MAw = '-' + MAw
        if B[2] == '0':
            MAu = u[k] + MAu
            k -= 1
        elif B[2] == '1':
            MAu = '-' + MAu
        B = backtrack[i+1, j+1, k+1]
    return int(S[-1,-1,-1]), MAv, MAw, MAu

def GreedySorting(P):
    approxReversalDistance = 0
    Plist = []
    for k in range(1, len(P)+1):
        if P[k-1] != k:
            index = [abs(i) for i in P].index(k)
            P = P[:k-1] + [-i for i in P][k-1: index+1][::-1] + P[index+1:]
            approxReversalDistance += 1
            Plist.append(AddSign(P))
        if P[k-1] == -k:
            P[k-1] = k
            approxReversalDistance += 1
            Plist.append(AddSign(P))
    return Plist

def AddSign(P):
    Psign = ''
    for p in P:
        if p > 0:
            Psign += '+' + str(p) + ' '
        elif p < 0:
            Psign += str(p) + ' '
    Psign = Psign[:-1]
    return Psign

def RemoveSign(P): ## i.e. P = '(+1 -2 -3)'
    P = P[1:-1]
    P = P.split(' ')
    Plist = []
    for p in P:
        if p[0] == '+':
            Plist.append(int(p[1:]))
        elif p[0] == '-':
            Plist.append(int(p))
    return Plist

def Breakpoints(P):
    breakpoints = 0
    P = [0] + P + [len(P) + 1]
    for i in range(len(P) - 1):
        if P[i+1] - P[i] != 1:
            breakpoints += 1
    return breakpoints

def ChromosomeToCycle(Chromosome): ## i.e. Chromosome = '(+1 -2 -3)'
    Nodes = []
    Chromosome = RemoveSign(Chromosome)
    for j in range(len(Chromosome)):
        i = Chromosome[j]
        if i > 0:
            Nodes.append(2*i - 1)
            Nodes.append(2*i)
        elif i < 0:
            Nodes.append(-2*i)
            Nodes.append(-2*i - 1)
    return Nodes

def ChromosomesToCycles(Chromosomes): ## i.e. Chromosomes = '(+1 -3 -6 -5)(+2 -4)'
    Cycles = []
    while len(Chromosomes) != 0:
        s = Chromosomes.index('(')
        e = Chromosomes.index(')')
        C = Chromosomes[s:e+1]
        cycle = ChromosomeToCycle(C)
        Cycles.append(cycle)
        Chromosomes = Chromosomes[e+1:]
    return Cycles

def CycleToChromosome(Nodes): ## i.e. Nodes = '[1, 2, 4, 3, 6, 5]'
    Chromosome = []
    for i in range(int(len(Nodes) / 2)):
        if Nodes[2*i] < Nodes[2*i + 1]:
            Chromosome.append(int(Nodes[2*i + 1] / 2))
        else:
            Chromosome.append(- int(Nodes[2*i] / 2))
    Chromosome = AddSign(Chromosome)
    return Chromosome

def ColoredEdges(P): ## i.e. P = '(+1 -2 -3)(+4 +5 -6)'
    Edges = set()
    cycles = ChromosomesToCycles(P)
    for cycle in cycles:
        for i in range(int(len(cycle) / 2) - 1):
            Edges.add((cycle[2*i + 1], cycle[2*i + 2]))
        Edges.add((cycle[-1], cycle[0]))
    Edges = [tuple(sorted(i)) for i in Edges]
    Edges = sorted(Edges)
    Edges = ', '.join([str(i) for i in Edges])
    return Edges

def Pair(i):
    if i % 2 == 0:
        return i - 1
    elif i % 2 == 1:
        return i + 1

def GraphToGenome(Edges): ## i.3. Edges = '(2, 4), (3, 6), (5, 1), (7, 9), (10, 12), (11, 8)'
    P = []
    Edges = Edges[1:-1]
    Edges = Edges.split('), (')
    Edges = [int(j) for i in Edges for j in i.split(', ')]
    #while Edges: ## Edges = [k2, k3, k4, ..., k1, k10, k11 ...]
    while Edges: ## Edges = [k0, k1, k2, ..., k2n]
        k0 = Edges[0]
        k1 = Edges[1]
        start = Pair(k0)
        cycle = [start]
        while k1 != start:
            cycle += [k0, k1]
            k0 = Pair(k1)
            ik0 = Edges.index(k0)
            if  ik0 % 2 == 0:
                ik1 = ik0 + 1
            elif ik0 % 2 == 1:
                ik1 = ik0 - 1
            k1 = Edges[ik1]
        else:
            cycle.append(k0)
            P.append(CycleToChromosome(cycle))
            Edges = [i for i in Edges if i not in cycle]
    P = ')('.join(P)
    P = '(' + P + ')'
    return P

def twoBreakDistance(P1, P2): ## i.e. P1 = '(+1 +2 +3 +4 +5 +6)'  P2 = '(+1 -3 -6 -5)(+2 -4)'
    redEdges = ColoredEdges(P1)
    blueEdges = ColoredEdges(P2)
    nBlocks = len(redEdges)
    Cycles = []
    while redEdges:
        re = list(redEdges)[0]
        start = re[0]
        be1 = re[0]         ## i.e. re = (re0, re1)  be = (be0, be1)
        cycle = [start]
        while True:
            re, re1 = Walk(redEdges, be1)
            cycle.append(re1)
            redEdges.remove(re)
            be, be1 = Walk(blueEdges, re1)
            if be1 == start:
                Cycles.append(cycle)
                break
            cycle.append(be1)
            blueEdges.remove(be)
    nCycles = len(Cycles)
    twobreakdistance = nBlocks - nCycles
    return twobreakdistance


def Walk(Edges, node_from):
    for edge in Edges:
        if edge[0] == node_from:
            node_to = edge[1]
            walked_edge = edge
        elif edge[1] == node_from:
            node_to = edge[0]
            walked_edge = edge
    return walked_edge, node_to

def twoBreakOnGenomeGraph(Edges, ilist): ## i.e. '(2, 4), (3, 8), (7, 5), (6, 1)'
                                                  ##      '1, 6, 3, 8'
    Edges = Edges[1:-1]
    Edges = Edges.split('), (')
    Edges = [int(j) for i in Edges for j in i.split(', ')]
    ilist = [int(i) for i in ilist.split(', ')]
    Edges = [i for i in Edges if i not in ilist]
    Edges += [ilist[i] for i in [0, 2, 1, 3]]
    Edges_set = set()
    for i in range(int(len(Edges) / 2)):
        Edges_set.add((Edges[2*i], Edges[2*i + 1]))
    Edges_set = [tuple(sorted(i)) for i in Edges_set]
    Edges_set = sorted(Edges_set)
    Edges = str(Edges_set)
    Edges = Edges[1:-1]
    return Edges

def twoBreakOnGenome(Edges, ilist): ## i.e. '(+1 -2 -4 +3)'
                                    ##      '1, 6, 3, 8'
    Edges = ColoredEdges(Edges) 
    Edges = twoBreakOnGenomeGraph(Edges, ilist) 
    Edges = GraphToGenome(Edges)
    return Edges

def StrToTupleSet(S):
    S = S[1:-1]
    S = S.split('), (')
    A = set()
    for s in S:
        l = tuple([int(i) for i in s.split(', ')])
        A.add(l)
    return A

def ShortestRearrangementScenario(P, Q):
    print(P)
    RedEdges = ColoredEdges(P)
    BlueEdges = ColoredEdges(Q)
    while RedEdges != BlueEdges: ## need to be sorted before compared
        R = StrToTupleSet(RedEdges)
        B = StrToTupleSet(BlueEdges)
        for be in B:
            if be not in R:
                i1, i2 = [be[0], Walk(R, be[0])[1]]
                i3, i4 = [be[1], Walk(R, be[1])[1]]
                ilist = [i1, i2, i3, i4]
                ilist = ', '.join([str(i) for i in ilist])
                RedEdges = twoBreakOnGenomeGraph(RedEdges, ilist)
                P = twoBreakOnGenome(P, ilist)
                print(P)
                break

                                                ##input
                                                ##'(+1 -2 -3 +4)'
                                                ##'(+1 +2 -4 -3)'

                                                ##output
                                                ##'(+1 -2 -3 +4)'
                                                ##'(+1 -2 -3)(+4)'
                                                ##'(+1 -2 -4 -3)'
                                                ##'(-3 +1 +2 -4)'
                                                
#def Sharedkmers(P, Q, k):
#    from Replication import ReverseComplement as RC
#    Shared = []
#    for i in range(len(P) - k + 1):
#        p = P[i:i+k]
#        for j in range(len(Q) - k + 1):
#            q = Q[j:j+k]
#            if p == q:
#                print((i, j))
#                Shared.append((i, j))
#            elif p == RC(q):
#                print((i, j))
#                Shared.append((i, j))
#    return Shared

def Sharedkmers(P, Q, k):
    from Replication import ReverseComplement as RC
    import sys
    sys.path.append('C:\\Users\\kimur\\Documents\\src\\Python_for_Genomic_Data_Science')
    from kmer_index import SubseqIndex as SI
    SIp = SI(P, k, 1)
    SIq = SI(Q, k, 1)
    kmers = set()
    Shared = []
    for i in range(len(SIp.index)):
        kmer = SIp.index[i][0]
        kmers.add(kmer)
    for kmer in kmers:
        Phits = SIp.query(kmer)
        Qhits = set(SIq.query(kmer) + SIq.query(RC(kmer)))
        Shared += list(itertools.product(Phits, Qhits))
        Shared.sort()
    return Shared



