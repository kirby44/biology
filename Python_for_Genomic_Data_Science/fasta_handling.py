#!/usr/bin/env python3

#fasta_path = "./Documents/dna.example.fasta"
#fasta_path = "./Documents/dna2.fasta"
def read_fasta(fasta_path):
    title_list = []
    sequence_list = []
    sequence = ""
    with open(fasta_path,'r') as fh:
        for line in fh:
            if line[0] == ">":
                title_list.append(line)
                sequence_list.append(sequence)
                sequence = ""
            else:
                sequence = sequence + line.rstrip()
        else:
            sequence_list.append(sequence)
    sequence_list.remove('')
    length_list = [len(s) for s in sequence_list] 
    max_length = max(length_list)
    max_title = title_list[length_list.index(max_length)]
    min_length = min(length_list)
    min_title = title_list[length_list.index(min_length)]
    return sequence_list,title_list

def readFastq(fastq_path):
    sequences = []
    qualities = []
    with open(fastq_path,'r') as fh:
        while True:
            fh.readline()
            seq = fh.readline().rstrip()
            fh.readline()
            qual = fh.readline().rstrip()
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences,qualities


#open reading frame(ORF)
#from textwrap import wrap
#[sequence_list,title_list] = read_fasta(fasta_path)
#start_codon = {"ATG"}
#stop_codon = {"TAA","TAG","TGA"}
def longest_orf_in_reading_frame(sequence,reading_frame):
    sequence = sequence[reading_frame - 1 :]
    codon_list = wrap(sequence,3)
    if len(codon_list) != 3:
        codon_list = codon_list[:-1]
    has_start = any([x in start_codon for x in codon_list])
    has_stop = any([x in stop_codon for x in codon_list])
    if not has_start or not has_stop:
        return 0,0,[]
    orf_list = []
    orf = []
    is_orf = False
    for i in range(len(codon_list)):
        codon = codon_list[i]
        if codon in start_codon:
            orf = [codon]
            is_orf = True
            continue
        if codon in stop_codon:
            orf.append(codon)
            orf_list.append(orf)
            orf = []
            is_orf = False
            continue
        if is_orf:
            orf.append(codon)
    length_list = [len(s) for s in orf_list]
    longest_length = max(length_list)
    longest_index = length_list.index(longest_length)
    longest_orf = orf_list[longest_index]
    longest_sequence = ""
    for codon in longest_orf:
        longest_sequence = longest_sequence + codon
    starting_point = sequence.index(longest_sequence)
    return longest_length,starting_point,longest_orf
def longest_orf_in_sequence(sequence):
    length_list = []
    starting_point_list = []
    orf_list = []
    for reading_frame in [1,2,3]:
    #for reading_frame in [3]:
        [longest_length,starting_point,longest_orf] = longest_orf_in_reading_frame(sequence,reading_frame)
        length_list.append(longest_length)
        starting_point_list.append(starting_point)
        orf_list.append(longest_orf)
    longest_length_in_sequence = max(length_list)
    longest_index = length_list.index(longest_length_in_sequence)
    longest_starting_point = starting_point_list[longest_index]
    longest_orf_in_sequence = orf_list[longest_index]
    return longest_length_in_sequence,longest_starting_point,longest_orf_in_sequence
def longest_orf_in_sequence_list(sequence_list):
    length_list = []
    starting_point_list = []
    orf_list = []
    for i in range(len(sequence_list)):
        [length,starting_point,orf] = longest_orf_in_sequence(sequence_list[i])
        length_list.append(length)
        starting_point_list.append(starting_point)
        orf_list.append(orf)
#        print(title_list[i])
#        print(length,starting_point,orf)
    longest_length = max(length_list)
    longest_index = length_list.index(longest_length)
    longest_starting_point = starting_point_list[longest_index]
    longest_orf = orf_list[longest_index]
#    print("longest length in this sequence list is ",longest_length)
#    print("starting point is ",longest_starting_point)
#    print("title is ",title_list[longest_index])
#    print("open reading frame is ",longest_orf)
#longest_orf_in_sequence_list(sequence_list)
