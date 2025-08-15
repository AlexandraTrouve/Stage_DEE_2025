from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("genes_files")
parser.add_argument("all_sequences")
parser.add_argument("filtered_sequences")
args = parser.parse_args()

# Read gene_ID from the txt file
with open(args.genes_files, "r") as f:
    ids_to_keep = set(line.strip() for line in f if line.strip())

# keep only sequences present in the previous txt file
written_ids = set()

with open(args.filtered_sequences, "w") as out_handle:
    for record in SeqIO.parse(args.all_sequences, "fasta"):
        full_gene_id = record.id
        gene_id = record.id.split(":")[-1]
        if gene_id in ids_to_keep and gene_id not in written_ids:
            record.description = ""
            SeqIO.write(record, out_handle, "fasta")
            written_ids.add(gene_id)
