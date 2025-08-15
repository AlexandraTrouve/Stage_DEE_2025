import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("evolution_scores")
parser.add_argument("mean_evolution_scores")
args = parser.parse_args()

input_file = args.evolution_scores
output_file = args.mean_evolution_scores


scores_dict = {}
with open(input_file, 'r') as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        parts = line.split("\t")
        id = parts[0]
        scores = [float(score) for score in parts[1:]]
        scores_dict[id] = scores

average_scores = {
    gene_id: (sum(scores) / len(scores) if scores else 0.0)
    for gene_id, scores in scores_dict.items()}

df = pd.DataFrame.from_dict(average_scores, orient='index', columns=['Score'])
df.index.name = 'Gene_ID'
df.to_csv(output_file)