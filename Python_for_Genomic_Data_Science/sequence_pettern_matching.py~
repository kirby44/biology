#!/usr/bin/env python3

#fasta_path = './Documents/fast/chr1.GRCh38.excerpt.fasta'
#fasta_path = '.\Documents\lambda_virus.fa'
fastq_path = './Documents/fast/ERR266411_1.for_asm.fastq'

from fasta_handling import read_fasta, readFastq
#sequence_list, title_list = read_fasta(fasta_path)
sequence_list, qualities = readFastq(fastq_path)

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

def overlap_all_pairs(reads, min_length=3):
    kmer_dic = {}
    for read in reads:
        for i in range(len(read) - min_length +1):
            kmer = read[i:i + min_length]
            if kmer in kmer_dic:
                kmer_dic[kmer].add(read)
            else:
                kmer_dic[kmer] = set([read])
    print('kmer_dic has been made')
    overlap_pairs = []
    for read in reads:
        suffix = read[-min_length:]
        b_candidates = kmer_dic[suffix]
        for b in b_candidates:
            if b != read:
                if overlap(read, b, min_length) > 0:
                    overlap_pairs.append([read,b])
    return overlap_pairs
#reads = ['ABCDEFG', 'EFGHIJ', 'HIJABC']
#reads = ['CGTACG', 'TACGTA', 'GTACGT', 'ACGTAC', 'GTACGA', 'TACGAT']
reads = sequence_list
print(overlap_all_pairs(reads, 30))

#p = 'GCGTATGC'
#t = 'TATTGGCTATACGGTT'
#p = 'GCTGATCGATCGTACG'
#p = 'GATTTACCAGATTGAG'
#t = sequence_list[0]
def editDistance(x, y):
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
    # Edit distance is the value in the bottom right corner of the matrix
    return D[-1][-1]

def approximate_editDistance(x, y):
    # Create distance matrix
    D = []
    for i in range(len(x)+1):
        D.append([0]*(len(y)+1))
    # Initialize first row and column of matrix
    for i in range(len(x)+1):
        D[i][0] = i
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
    return min(D[-1])
#print(approximate_editDistance(p, t))

def reverseComplement(s):
    complement = {'A':'T','C':'G','G':'C','T':'A','N':'N'}
    rc = ''
    for chr in s:
        rc = complement[chr] + rc
    return rc

def naive_with_rc(p,t):
    rc = reverseComplement(p)
    occurrences = []
    for i in range(len(t)-len(p)+1):
        match_p = True
        match_rc = True
        for j in range(len(p)):
            if p[j] != t[i+j]:
                match_p = False
                break
        for j in range(len(p)):
            if rc[j] != t[i+j]:
                match_rc = False
                break
        if match_p or match_rc:
            occurrences.append(i)
    return occurrences
#o = naive_with_rc(p,t)

def naive_with_counts(p,t):
    occurrences = []
    num_comparison = 0
    num_alignment = len(t) - len(p) + 1
    for i in range(len(t)-len(p)+1):
        match_p = True
        for j in range(len(p)):
            num_comparison += 1
            if p[j] != t[i+j]:
                match_p = False
                break
        if match_p:
            occurrences.append(i)
    return occurrences, num_alignment, num_comparison

def naive_with_mm(p,t):
    occurrences = []
    for i in range(len(t)-len(p)+1):
        match = True
        n_mismatch = 0
        for j in range(len(p)):
            if p[j] != t[i+j]:
                n_mismatch += 1
                if n_mismatch > 2:
                    match = False
                    break
        if match:
            occurrences.append(i)
    return occurrences
#o = naive_with_mm(p,t)

def boyer_moore(p, p_bm, t):
    occurrences = []
    i = 0
    while i < len(t) - len(p):
        shift = 1
        mismatched = False
        for j in range(len(p)-1,-1,-1):
            if p[j] != t[i+j]:
                s_bc = p_bm.bad_character_rule(j,t[i+j])
                s_gs = p_bm.good_suffix_rule(j)
                shift = max(s_bc,s_gs,shift)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            s_gs = p_bm.match_skip()
            shift = max(s_gs,shift)
        i += shift
    return occurrences

def boyer_moore_with_counts(p, p_bm, t):
    occurrences = []
    i = 0
    num_comparison = 0
    num_alignment = 0
    while i < len(t) - len(p) + 1:
        num_alignment += 1
        shift = 1
        mismatched = False
        for j in range(len(p)-1,-1,-1):
            num_comparison += 1
            if p[j] != t[i+j]:
                s_bc = p_bm.bad_character_rule(j,t[i+j])
                s_gs = p_bm.good_suffix_rule(j)
                shift = max(s_bc,s_gs,shift)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            s_gs = p_bm.match_skip()
            shift = max(s_gs,shift)
        i += shift
    return occurrences, num_alignment, num_comparison

def approximate_match(p, t, n):
    import kmer_index
    index = kmer_index.Index(t,8)
    segment_length = round(len(p) / (n+1))
    all_matches = set()
    total_hits = 0
    for i in range(n+1):
        start = i*segment_length
        end = min((i+1)*segment_length,len(p))
        segment_p = p[start:end]
        matches = index.query(segment_p)
        total_hits += len(matches)
        print(len(matches),segment_p)
        for m in matches:
            if m < start or m-start+len(p) > len(t):
                continue
            mismatches = 0
            for j in range(0,start):
                if p[j] != t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            for j in range(end,len(p)):
                if p[j] != t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            if mismatches <= n:
                all_matches.add(m-start)
    print(total_hits)
    return list(all_matches)

#import bm_preproc
#from kmer_index import SubseqIndex
##t = 'to-morrow and to-morrow and to-morrow creeps in this petty pace'
##p = 'to-morrow and to-morrow '
#t = sequence_list[0]
#p  = 'GGCGCGGTGGCTCACGCCTGTAAT'
#subseq_ind = SubseqIndex(t, 8, 3)
#occurrences, num_index_hits = subseq_ind.query_subseq(p)
#print(occurrences)
#print(num_index_hits)

##p  = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
#p  = 'GGCGCGGTGGCTCACGCCTGTAAT'
#p1 = 'GGCGCGGT'
#p2 = 'GGCTCACG'
#p3 = 'CCTGTAAT'
#t = sequence_list[0]
#matches = approximate_match(p, t, 2)
#print(len(matches), matches)
#
#p_bm = bm_preproc.BoyerMoore(p1)
#occurrences, num_alignments, num_character_comparisons = boyer_moore_with_counts(p1, p_bm, t)
#print(len(occurrences), occurrences, num_alignments, num_character_comparisons)
#p_bm = bm_preproc.BoyerMoore(p2)
#occurrences, num_alignments, num_character_comparisons = boyer_moore_with_counts(p2, p_bm, t)
#print(len(occurrences), occurrences, num_alignments, num_character_comparisons)
#p_bm = bm_preproc.BoyerMoore(p3)
#occurrences, num_alignments, num_character_comparisons = boyer_moore_with_counts(p3, p_bm, t)
#print(len(occurrences), occurrences, num_alignments, num_character_comparisons)

#p = 'word'
#t = 'there would have been a time for such a word'
#alphabet = 'abcdefghijklmnopqrstuvwxyz '
#occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, t)
#print(occurrences, num_alignments, num_character_comparisons)
#p_bm = bm_preproc.BoyerMoore(p, alphabet)
#occurrences, num_alignments, num_character_comparisons = boyer_moore_with_counts(p, p_bm, t)
#print(occurrences, num_alignments, num_character_comparisons)
#
#p = 'needle'
#t = 'needle need noodle needle'
#occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, t)
#print(occurrences, num_alignments, num_character_comparisons)
#p_bm = bm_preproc.BoyerMoore(p, alphabet)
#occurrences, num_alignments, num_character_comparisons = boyer_moore_with_counts(p, p_bm, t)
#print(occurrences, num_alignments, num_character_comparisons)
#
