import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb
import argparse

cluster = argparse.ArgumentParser()
cluster.add_argument("pred")
cluster.add_argument("scores")
args = cluster.parse_args()


input_f = pd.read_csv(args.pred)
input_f = input_f.drop(columns = ['Kmer'])

# correlation matrix construction
correlation_matrix = input_f.corr()
correlation_matrix

# heatmap construction
plt.figure(figsize=(6, 5))
sb.heatmap(correlation_matrix,
            annot=True,  
            fmt=".2f",
            linewidths=.5,
            cbar_kws={'label': 'Correlation coefficient'})

plt.title('Correlation Heatmap between Gene Clusters')
plt.xlabel('Gene clusters')
plt.ylabel('Gene clusters') 
plt.tight_layout()
plt.savefig(args.scores)
plt.show()
plt.close()