##!/user/bin/env python3

#fasta_path = "./Documents/dna.example.fasta"
fasta_path = "./Documents/dna2.fasta"
import sys
sys.path.append("\\Users\\kimur\\Documents\\src")
import fasta_handling
[sequence_list,title_list] = fasta_handling.read_fasta(fasta_path)
#fasta_handling.longest_orf_in_sequence_list(sequence_list)

def repeat_per_length(sequence_list,repeat_length):
    repeat_per_str_dic = {}
    for sequence in sequence_list:
        for i in range(len(sequence)-repeat_length+1):
            str_i = sequence[i:i + repeat_length]
            if str_i in repeat_per_str_dic:
                repeat_per_str_dic[str_i] += 1
            else:
                repeat_per_str_dic[str_i] = 1
    return repeat_per_str_dic
repeat_dic_list = [repeat_per_length(sequence_list,i+1) for i in range(15)]

#6
dic6 = repeat_dic_list[5]
max_key = max(dic6,key=dic6.get)
max_r = dic6[max_key]
print(max_r)

#12
dic12 = repeat_dic_list[11]
max_key = max(dic12,key=dic12.get)
max_r = dic12[max_key]
count = 0
for key in dic12.keys():
    if dic12[key] == max_r:
        print(key)
        count += 1
print(count)

#7
dic7 = repeat_dic_list[6]
max_key = max(dic7,key=dic7.get)
max_r = dic7[max_key]
print(max_key)

print("finish")

