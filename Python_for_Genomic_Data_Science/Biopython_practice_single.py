#!/usr/bin/env python3

from Bio.Blast import NCBIWWW,NCBIXML

fasta_path = "./Documents/dna.example.single.fasta"
fasta_string = open(fasta_path).read()
result_handle = NCBIWWW.qblast("blastn","nt",fasta_string)
blast = NCBIXML.read(result_handle)

for alignment in blast.alignments:
    for hsp in alignment.hsps:
        print(alignment.title)
        print(alignment.length)
        print(hsp.expect)
        print(hsp.query)
        print(hsp.match)
        print(hsp.sbjct)
