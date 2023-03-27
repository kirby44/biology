#!/usr/bin/env python3

fastq_path = './Documents/fast/ads1_week4_reads.fq'

from fasta_handling import readFastq
sequence_list, qualities = readFastq(fastq_path)

def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's suffx in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match

import itertools

def scs(ss):
    """ Returns shortest common superstring of given
        strings, which must be the same length """
    shortest_sup = None
    len_ss = len(''.join(ss))
    k = range(1,len_ss+1)
    v = [0] * len_ss
    num_sup_dic = dict(zip(k,v))
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss)-1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            # add non-overlapping portion of B to superstring
            sup += ssperm[i+1][olen:]
        num_sup_dic[len(sup)] += 1
        if shortest_sup is None or len(sup) < len(shortest_sup):
            shortest_sup = sup  # found shorter superstring
    return shortest_sup, num_sup_dic  # return shortest

#ss = ['CCT', 'CTT', 'TGC', 'TGG', 'GAT', 'ATT']
#SCS, num_sup_dic = scs(ss)
#print(SCS)
#print(len(SCS))
#print(num_sup_dic)

def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's prefix in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match

def overlap_graph(reads):
    nodes = set(reads)
    edges = []
    current_max_overlap = 0
    for pair in itertools.permutations(nodes,2):
        o = overlap(pair[0],pair[1])
        if o > current_max_overlap:
            current_max_overlap = o
            edges = [pair]
        elif o == current_max_overlap:
            edges.append(pair)
    return nodes, edges, current_max_overlap

def greedy_scs(reads):
    nodes, edges, current_max_overlap = overlap_graph(sequence_list)
    #print("overlap graph has been created")
    while len(edges) > 0:
        a,b = edges[0]
        nodes.remove(a)
        nodes.remove(b)
        new_node = a + b[current_max_overlap:]
        if len(nodes) == 0:
            return new_node
        del_edges = set()
        for edge in edges:
            if a in edge or b in edge:
                del_edges.add(edge)
        for edge in del_edges:
            edges.remove(edge)
        for node in nodes:
            new_overlap_prefix = overlap(new_node,node)
            assert new_overlap_prefix <= current_max_overlap, 'new nodes excced current max overlap'
            if new_overlap_prefix == current_max_overlap:
                edges.append((new_node, node))
            new_overlap_suffix = overlap(node,new_node)
            assert new_overlap_suffix <= current_max_overlap, 'new nodes excced current max overlap'
            if new_overlap_suffix == current_max_overlap:
                edges.append((node, new_node))
        nodes.add(new_node)
        #print('num nodes = ', len(nodes),'overlap = ', current_max_overlap)
        if len(edges) == 0:
            #print('creating overlap graph')
            nodes, edges, current_max_overlap = overlap_graph(nodes)
scs = greedy_scs(sequence_list)
#print(len(scs))
print(scs)
