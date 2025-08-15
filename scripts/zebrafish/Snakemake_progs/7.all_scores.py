from Bio import SeqIO
import pandas as pd
import argparse
from Bio.Seq import Seq


# Get SVM score for each kmer of a given model
def get_svm_dict(path_model):
    svm_dict = {}
    with open(path_model, 'r') as model:
        for i in model.readlines():
            i = i.strip("\n")
            i = i.split("\t")
            kmer = i[0]
            svm_score = float(i[1])
            rev_kmer = str(Seq(kmer).reverse_complement())  # add the reverse complement kmer
            svm_dict[kmer] = svm_score
            svm_dict[rev_kmer] = svm_score
    return svm_dict

def calculate_svm(seq, svm_dict, kmer_len=10):
    svm = 0
    for pos in range(len(seq) - kmer_len + 1):  # sliding window of kmer length
        kmer = seq[pos:pos + kmer_len]
        if kmer in svm_dict:
            svm += svm_dict[kmer]  # sum of SVM for each kmer
        else:
            continue
    return round(svm, 7)

parser = argparse.ArgumentParser()
parser.add_argument("scores")
parser.add_argument("Sequences")
parser.add_argument("sequences_scores")
args = parser.parse_args()

fasta_file = args.Sequences
scores_file = args.scores
output_file = args.sequences_scores
all_svm = get_svm_dict(scores_file)

all_scores = {}

for record in SeqIO.parse(fasta_file, "fasta"):
    gene_id = record.id
    sequence = str(record.seq)
    Sequence = sequence.upper()
    score = calculate_svm(Sequence, all_svm, kmer_len=10)
    all_scores[gene_id] = score

df = pd.DataFrame.from_dict(all_scores, orient='index', columns=['SVM_Score'])
df.index.name = 'Gene_ID'
df.to_csv(output_file)