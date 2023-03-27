import numpy as np
import itertools
import math
import collections
from copy import deepcopy

def StringsToEdges(path): ## path = file path to strings of edges
    Edges = dict()
    with open(path, 'r') as f:
        lines = f.read().splitlines()[1:]                       ## assume head line is n_leaves
        for line in lines:
            key, value = line.split(':')
            value = int(value)
            Edges[key] = value
    return Edges

def Path(Edges, start, end, innode = ''):
    start = str(start)
    end = str(end)
    outnodes = [i.split('->')[1] for i in Edges.keys() if i.split('->')[0] == start]
    if innode in outnodes:  ## exclude return path  i.e. 2 -> 3 -> 2
        outnodes.remove(innode)
    if end in outnodes:
        return [start, end]
    else:
        for o in outnodes:
            p = Path(Edges, o, end, start)
            if p != None:
                if p[-1] == end:
                    return [start] + p

def Distance(Edges, a, b):
    a = str(a)
    b = str(b)
    if a == b:
        return 0
    d = 0
    p = Path(Edges, a, b, '')
    for i in range(len(p) - 1):
        key = p[i] + '->' + p[i+1]
        d += Edges[key]
    return d

def Matrix(Edges_path, n_leaves):
    M = np.zeros((n_leaves, n_leaves))
    Edges = StringsToEdges(Edges_path)
    for a, b in itertools.combinations(range(n_leaves), 2):
        d = Distance(Edges, str(a), str(b))
        M[a, b] = d
    M = M + np.transpose(M)
    return M


                                    ##      0->4:11
                                    ##      1->4:2
                                    ##      2->5:6
                                    ##      3->5:7
                                    ##      4->0:11
                                    ##      4->1:2
                                    ##      4->5:4
                                    ##      5->4:4
                                    ##      5->3:7
                                    ##      5->2:6
                                    
                                    ##  0	13	21	22
                                    ##  13	0	12	13
                                    ##  21	12	0	13
                                    ##  22	13	13	0
def StrToMatrix(path):
    M = np.loadtxt(path)

def Limb(D, j):
    Limb = math.inf
    n = len(D)
    l = [i for i in range(n) if i != j]
    for i, k in itertools.combinations(l, 2):
        dij = D[i, j]
        djk = D[j, k]
        dik = D[i, k]
        d = (dij + djk - dik) / 2
        if Limb > d:
            Limb = d
    return Limb

def Increment(T, n):
    IncrementedT = dict()
    for key, value in T.items():
        innode, outnode = [int(i) for i in key.split('->')]
        if innode not in range(n-1):
            innode += 1
        if outnode not in range(n-1):
            outnode += 1
        key = str(innode) + '->' + str(outnode)
        IncrementedT[key] = value
    return IncrementedT

def nNodes(T):
    keys = T.keys()
    nodes = [k.split('->') for k in keys]
    innodes = [n[0] for n in nodes]
    outnodes = [n[1] for n in nodes]
    nNodes = len(set(innodes + outnodes))
    return nNodes

def Sort(T):
    SortedT = dict()
    n_list = []
    for k, v in T.items():
        nodes = tuple([int(i) for i in k.split('->')] + [v])
        n_list.append(nodes)
    n_list.sort()
    for n in n_list:
        innode = str(n[0])
        outnode = str(n[1])
        if isinstance(n[2], int):
            v = int(n[2])
        elif isinstance(n[2], float):
            v = format(n[2], '.3f')
        key = innode + '->' + outnode
        SortedT[key] = v
    return SortedT

def AdditivePhylogeny(D):
    n = len(D)
    if n == 2:
        return {'0->1':D[0, 1], '1->0':D[1, 0]}
    LimbLength = Limb(D, n-1)
    D[n-1, 0:-1] -= LimbLength
    D[0:-1, n-1] -= LimbLength
    for i, k in itertools.combinations(range(n-1), 2):
        if D[i, k] == D[i, n-1] + D[n-1, k]:
            x = D[i, n-1]
            break
    D = D[:-1,:-1]
    T = AdditivePhylogeny(D)
    if n > 3:
        T = Increment(T, n)

    ## find edge e which has attatchment point v
    p = Path(T, i, k)
    for l in range(1, len(p)):
        e0 = p[l-1]
        e1 = p[l]
        die0 = Distance(T, i, e0)
        die1 = Distance(T, i, e1)
        if x < die1:
            ## divide an old edge at v into two new edges
            v = str(nNodes(T) + 1)
            T[e0 + '->' + v] = x - die0
            T[v + '->' + e0] = x - die0
            T[e1 + '->' + v] = T[e0 + '->' + e1] - T[e0 + '->' + v]
            T[v + '->' + e1] = T[e0 + '->' + e1] - T[e0 + '->' + v]
            del T[e0 + '->' + e1]
            del T[e1 + '->' + e0]
            break
        elif x == die1:
            v = e1
            break
    
    ## add edge v -> n
    T[v + '->' + str(n-1)] = LimbLength
    T[str(n-1) + '->' + v] = LimbLength
    return T

def UPGMA(D):
    T = dict()
    Clusters = list(range(len(D)))
    Age = {k:0 for k in range(len(D))}
    D[np.tril_indices(len(D))] = np.inf
    while len(D) > 1:
        w = np.where(D == np.min(D))
        i = w[0][0]
        j = w[1][0]
        Cnew = len(Clusters) 
        Age[Cnew] = D[i, j] / 2
        T[str(Cnew) + '->' + str(Clusters[i])] = Age[Cnew] - Age[Clusters[i]]
        T[str(Clusters[i]) + '->' + str(Cnew)] = Age[Cnew] - Age[Clusters[i]]
        T[str(Cnew) + '->' + str(Clusters[j])] = Age[Cnew] - Age[Clusters[j]]
        T[str(Clusters[j]) + '->' + str(Cnew)] = Age[Cnew] - Age[Clusters[j]]
        newcol = [] 
        for x in range(len(D)):
            if x not in (i, j):
                dix = min(D[i, x], D[x, i])
                djx = min(D[j, x], D[x, j])
                dnewx = (dix + djx) / 2
                newcol.append(dnewx)
        D = np.delete(D, [i, j], 0)
        D = np.delete(D, [i, j], 1)
        D = np.insert(D, len(D), newcol, 1)
        D = np.insert(D, len(D), np.inf, 0)
        Clusters.remove(Clusters[j])
        Clusters.remove(Clusters[i])
        Clusters.append(Cnew)
    return T
    
def AddEdge(T, innode, outnode, value):
    key1 = str(innode) + '->' + str(outnode)
    key2 = str(outnode) + '->' + str(innode)
    T[key1] = value
    T[key2] = value
    return T

def NeighborJoining(D, Clusters = []):
    T = dict()
    n = len(D)
    if Clusters == []:
        Clusters = list(range(n))
    if n == 2:
        T = AddEdge(T, Clusters[0], Clusters[1], D[0, 1])
        return T
    D[np.tril_indices(n)] = np.inf
    Dnj = np.zeros((n, n))
    Dnj[np.tril_indices(n)] = np.inf
    TD = []
    for x in range(n):
        tdx = sum(D[x][np.where(D[x] != np.inf)]) + sum(D[np.where(D[:,x] != np.inf), x][0]) 
        TD.append(tdx)
    for i, j in itertools.combinations(range(n), r = 2):
        Dnj[i, j] = (n - 2) * D[i, j] - TD[i] - TD[j]
    w = np.where(Dnj == np.min(Dnj))
    i = w[0][0]
    j = w[1][0]
    Ci = Clusters[i]
    Cj = Clusters[j]
    delta = (TD[i] - TD[j]) / (n - 2)
    Li = (D[i, j] + delta) / 2
    Lj = (D[i, j] - delta) / 2
    newcol = [] 
    for x in range(n):
        if x not in (i, j):
            dix = min(D[i, x], D[x, i])
            djx = min(D[j, x], D[x, j])
            dnewx = (dix + djx - D[i, j]) / 2
            newcol.append(dnewx)
    D = np.delete(D, [i, j], 0)
    D = np.delete(D, [i, j], 1)
    D = np.insert(D, len(D), newcol, 1)
    D = np.insert(D, len(D), np.inf, 0)
    Cnew = max(Clusters) + 1
    Clusters.remove(Clusters[j])
    Clusters.remove(Clusters[i])
    Clusters.append(Cnew)
    T = NeighborJoining(D, Clusters)
    T = AddEdge(T, Cnew, Ci, Li)
    T = AddEdge(T, Cnew, Cj, Lj)
    return T

def SubAssignCharacter(son, Cparent, S, Character = ['A', 'C', 'G', 'T']):
    if S[son, Cparent] > min(S[son]) + 1:
        minC = np.where(S[son] == min(S[son]))[0][0]
        return Character[minC]
    else:
        return Character[Cparent]

def AssignCharacter(Edges, nodes, S, Character = ['A', 'C', 'G', 'T']):
    while '' in nodes.values():
        for v in nodes.keys():
            if v not in Edges.keys():
                continue
            else:
                daughter, son = Edges[v]
            if nodes[daughter] != '':
                continue
            Cparent = Character.index(nodes[v])
            nodes[daughter] = SubAssignCharacter(daughter, Cparent, S, Character)
            nodes[son] = SubAssignCharacter(son, Cparent, S, Character)
    return nodes

def SingleSmallParsimony(Edges, nodes, Character = ['A', 'C', 'G', 'T']):
    Tag = {k:0 if v == '' else 1 for k,v in nodes.items()}
    S = np.full((len(nodes.keys()), len(Character)), np.inf)
    for v, ripe in Tag.items():
        if ripe == 1:
            for c in range(len(Character)):
                if nodes[v] == Character[c]:
                    S[v][c] = 0
                else:
                    S[v][c] = math.inf
    while 0 in Tag.values():
        for v, ripe in Tag.items():
            if ripe == 0:
                daughter, son = Edges[v]
                if Tag[daughter] == 1 and Tag[son] == 1:
                    Tag[v] = 1
                    mind = min(S[daughter])
                    mindc = np.where(S[daughter] == mind)[0]
                    mins = min(S[son])
                    minsc = np.where(S[son] == mins)[0]
                    for c in range(len(Character)):
                        if c in mindc:
                            alphad = 0
                        else:
                            alphad = 1
                        if c in minsc:
                            alphas = 0
                        else:
                            alphas = 1
                        svc = mind + alphad + mins + alphas
                        S[v, c] = svc
    print(S)
    root = max(nodes.keys())
    Sroot = min(S[root])
    Croot = np.where(S[root] == Sroot)[0][0]
    nodes[root] = Character[Croot]
    nodes = AssignCharacter(Edges, nodes, S, Character)
    return Sroot, nodes

def AddEdge(Edges, innode, outnode):
    if innode in Edges.keys():
        Edges[innode].append(outnode)
    else:
        Edges[innode] = [outnode]
    return Edges

def SmallParsimony(T, Character = ['A', 'C', 'G', 'T']):
    Score = 0
    n = int(T[0])
    Edges = dict()
    nodes = dict()
    for i in range(n):
        innode, string = T[i + 1].split('->')
        innode = int(innode)
        Edges = AddEdge(Edges, innode, i)
        nodes[i] = string
    for i in range(len(T) - n - 1):
        Edge = T[n + 1 + i]
        innode, outnode = [int(i) for i in Edge.split('->')]
        Edges = AddEdge(Edges, innode, outnode)
        nodes[innode] = ''
        nodes[outnode] = ''
    print(Edges)
    print('len(nodes) = ', len(nodes))
    for i in range(len(string)):
        nodes_i = {k:v[i] if len(v) == len(string) else '' for k,v in nodes.items()}
        Score_i, nodes_i = SingleSmallParsimony(Edges, nodes_i)
        print('i = ', i, 'Score = ', Score_i, nodes_i)
        print('')
        Score += Score_i
        for v in range(len(nodes)):
            if not len(nodes[v]) == len(string):
                nodes[v] += nodes_i[v] 
    print('Total score = ', Score)
    print(nodes)
                                ##  4
                                ##  4->CAAATCCC
                                ##  4->ATTGCGAC
                                ##  5->CTGCGCTG
                                ##  5->ATGGACGA
                                ##  6->4
                                ##  6->5

