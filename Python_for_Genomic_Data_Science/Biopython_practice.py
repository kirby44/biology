#!/usr/bin/env python3

from Bio.Blast import NCBIWWW,NCBIXML

fasta_path = "./Documents/dna.example.fasta"
fasta_string = open(fasta_path).read()
print("There are ",fasta_string.count(">"),"records")
print("qblast start")
result_handle = NCBIWWW.qblast("blastn","nt",fasta_string)
print("qblast end")

blast_records = NCBIXML.parse(result_handle)
for blast_record in blast_records:
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            print(alignment.title)
