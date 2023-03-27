#!/usr/bin/env python3

fastq_path = './Documents/fast/ERR037900_1.first1000.fastq'
from fasta_handling import readFastq
sequences, qualities = readFastq(fastq_path)

def qualities_to_decimal(qualities):
    qualities_decimal_list = []
    for q in qualities:
        qualities_decimal = [ord(i) - 33 for i in q]
        qualities_decimal_list.append(qualities_decimal)
    return qualities_decimal_list

import matplotlib.pyplot as plt
h_list = qualities_to_decimal(qualities)
for h in h_list:
    min_q = min(h)
    print(h.index(min_q))
    plt.plot(range(len(h)),h)
    plt.show()

