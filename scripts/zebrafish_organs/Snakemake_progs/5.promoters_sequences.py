#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 21 14:53:00 2025

@author: children
"""

import Bio
from Bio import SeqIO
import csv
import argparse

### genome conversion from ensembl to UCSC notation
file = "/Users/children/Desktop/Stage_2025/data/zebrafish_organs/sequences/danRer11.fa"
genome = SeqIO.to_dict(SeqIO.parse(file, "fasta"))

### fatsa files contaning promoters sequences
cluster = argparse.ArgumentParser()
cluster.add_argument("Set")
cluster.add_argument("Sequences")
args = cluster.parse_args()


promoters = []
with open(args.Set, 'r') as bed_file:
    for line in bed_file:
        parts = line.strip().split('\t')
        if len(parts) == 4:
            chromosome = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            gene = parts[3]
            seq = genome[chromosome].seq      # recovery of the chromosome associated sequence
            promoter_sequence = seq[(start - 1):end]   # extraction of the promoter sequence according to defined coordinates
            record_promoter = SeqIO.SeqRecord(promoter_sequence, id=gene)   # storage of the recovery promoter sequence
            promoters.append(record_promoter)
        else:
            chromosome = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            seq = genome[chromosome].seq      # recovery of the chromosome associated sequence
            promoter_sequence = seq[(start - 1):end]   # extraction of the promoter sequence according to defined coordinates
            record_promoter = SeqIO.SeqRecord(promoter_sequence)   # storage of the recovery promoter sequence
            promoters.append(record_promoter)
    output_file_name = args.Sequences
    SeqIO.write(promoters, output_file_name, "fasta")   # creation of the fasta file containing promoters sequences





