##!/user/bin/env python3

#fasta_path = "./Documents/dna.example.fasta"
fasta_path = "./Documents/dna2.fasta"
import sys
sys.path.append("\\Users\\kimur\\Documents\\src")
import fasta_handling
[sequence_list,title_list] = fasta_handling.read_fasta(fasta_path)
fasta_handling.longest_orf_in_sequence_list(sequence_list)

def detect_repeat(sequence,repeat_length):
    repeat_sequences = set()
    for i in range(len(sequence) - repeat_length + 1):
        repeat_sequences.add(sequence[i:i + repeat_length])
    max_sequence = ""
    max_repeat_times = 0
    repeat_start_point = 0
    for repeat_sequence in repeat_sequences:
        indices = [i for i in range(len(sequence) - repeat_length + 1) if sequence[i:i + repeat_length] == repeat_sequence] 
        if len(indices) == 1:
            continue
        for overlap in range(repeat_length):
            [max_o,index] = detect_repeat_per_overlap(indices,repeat_length,overlap)
            if max_repeat_times < max_o:
                max_sequence = repeat_sequence
                max_repeat_times = max_o
                repeat_start_point = index
    return max_sequence,max_repeat_times,repeat_start_point

def detect_repeat_per_overlap(indices,repeat_length,overlap):
    max_repeat_times_per_overlap = 0
    repeat_start_point = 0
    mod = repeat_length - overlap
    for remainder in range(mod):
        [max_m,index] = detect_repeat_per_remainder(indices,repeat_length,overlap,remainder)
        if max_repeat_times_per_overlap < max_m:
            max_repeat_times_per_overlap = max_m
            repeat_start_point = index
    return max_repeat_times_per_overlap,repeat_start_point
def detect_repeat_per_remainder(indices,repeat_length,overlap,remainder):
    mod = repeat_length - overlap
    mod_group = [i for i in indices if i % mod == remainder]
    if len(mod_group) < 2:
        return 0,0
    repeat_times_list = []
    for i in range(len(mod_group) - 1):
        k = 1
        repeat_times = 1
        while mod_group[i] + repeat_length * k == mod_group[i+k] and i+k < len(mod_group):
            repeat_times += 1
            k += 1
            if i+k == len(mod_group):
                repeat_times_list.append(repeat_times)
                break
        else:
            repeat_times_list.append(repeat_times)
    if repeat_times_list == []:
        return 0,0
    max_repeat_times_per_remainder = max(repeat_times_list)
    index_in_mod_group = repeat_times_list.index(max_repeat_times_per_remainder)
    index = mod_group[index_in_mod_group]
    return max_repeat_times_per_remainder,index

max_t = ""
max_s = ""
max_r = 1
max_p = 0
for i,sequence in enumerate(sequence_list):
    print(title_list[i])
    max_sequence = ""
    max_repeat_times = 1
    max_start_point = 0
#    for k in range(len(sequence)):
    for k in range(20):
        [s,r,p] = detect_repeat(sequence,k+1)
        print(s,r,p)
        if max_repeat_times < r:
            max_sequence = s
            max_repeat_times = r
            max_start_point = p
        if len(sequence) / max_repeat_times < k:
            break
    print(max_sequence,max_repeat_times,max_start_point)
    if max_r < max_repeat_times:
        max_t = title_list[i]
        max_s = max_sequence
        max_r = max_repeat_times
        max_p = max_start_point
print("--Max--")
print(max_t)
print(max_s,max_r,max_p)
        
        
