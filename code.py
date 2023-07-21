import os, re
import glob
import pandas as pd
import numpy as np
from collections import Counter
from Bio import SeqIO                                                                 
from Bio import SeqIO
print("done")

fasta_file = "Original_seq.fasta" # Input fasta file
wanted_file = "ID.txt" # Input interesting sequence IDs
result_file = "extracted_seq.fasta" # extracted seq
remove_file = "rmf.txt"
fasta_file1 = "primery_removedSeq.fasta" # file dose not have extracted seq
inserted_seq = "insertedSeq.fasta"
proc_prime = "final_processed_primery.fasta"

### file Containing ID which is to be extracted
wanted = set()
with open(wanted_file) as f:
    for line in f:
        line = line.strip()
        if line != "":
            wanted.add(line)

fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
with open(result_file, "w") as f:
    for seq in fasta_sequences:
        if seq.id in wanted:
            SeqIO.write([seq], f, "fasta")

### removed extracted seq from primery fasta file
remove = set()
with open(remove_file) as f:
    for line in f:
        line = line.strip()
        if line != "":
            remove.add(line)

fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')

with open(fasta_file1, "w") as f:
    for seq in fasta_sequences:
        nam = seq.id
        nuc = str(seq.seq)
        if nam not in remove and len(nuc) > 0:
            SeqIO.write([seq], f, "fasta")


### convert original fasta file into BED file
with open('Original_seq.fasta') as in_f, open('sequences.bed','w') as out_f:
    for record in SeqIO.parse(in_f, 'fasta'):
        out_f.write('{}\t0\t{}\n'.format(record.id, len(record)))
        
        
### reverse complement of extracted sequence
from Bio.SeqRecord import SeqRecord

def make_rc_record(record):
    return SeqRecord(seq = record.seq.reverse_complement(), \
                 id = "rc_" + record.id, \
                 description = "reverse complement")
records = map(make_rc_record, SeqIO.parse("extracted_seq.fasta", "fasta"))
SeqIO.write(records, "revComp.fasta", "fasta")

### Introduce 1 base in extrected seq
fin = open("extracted_seq.fasta", "rt")

fout = open("processedSeq.fasta", "wt")
for line in fin:
    #read replace the string and write to output file
    fout.write(line.replace('T', 'A'))
fin.close()
fout.close()

### Insert processed Seq into primery seq 
with open(proc_prime, 'w') as output_fas:
    file_count = 0
    for f in os.listdir():
        if f.startswith("processedSeq") or f.startswith("primery_removedSeq"):
            file_count += 1
            with open(f) as fh: 
                output_fas.writelines(fh)
                
