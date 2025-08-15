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
file = "/Users/children/Desktop/Stage_2025/Data/Drosophila_melanogaster.BDGP6.46.dna_sm.toplevel.fa"
ensembl_genome = SeqIO.to_dict(SeqIO.parse(file, "fasta"))
matchings = "/Users/children/Desktop/Stage_2025/Data/chromosome_correspondence_Ensembl2UCSC.txt"

correspondances = {}
with open(matchings, 'r') as txt_file:
    for ligne in txt_file:
        parts = ligne.strip().split('\t')
        ensembl_name = parts[0]
        UCSC_name = parts[1]
        correspondances[ensembl_name] = UCSC_name
UCSC_genome = {}
for ensembl_name, sequence in ensembl_genome.items():
    UCSC_name = correspondances[ensembl_name]
    UCSC_genome[UCSC_name] = sequence


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
            seq = UCSC_genome[chromosome].seq      # recovery of the chromosome associated sequence
            promoter_sequence = seq[(start - 1):end]   # extraction of the promoter sequence according to defined coordinates
            record_promoter = SeqIO.SeqRecord(promoter_sequence, id=gene)   # storage of the recovery promoter sequence
            promoters.append(record_promoter)
        else:
            chromosome = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            seq = UCSC_genome[chromosome].seq      # recovery of the chromosome associated sequence
            promoter_sequence = seq[(start - 1):end]   # extraction of the promoter sequence according to defined coordinates
            record_promoter = SeqIO.SeqRecord(promoter_sequence)   # storage of the recovery promoter sequence
            promoters.append(record_promoter)
    output_file_name = args.Sequences
    SeqIO.write(promoters, output_file_name, "fasta")   # creation of the fasta file containing promoters sequences





