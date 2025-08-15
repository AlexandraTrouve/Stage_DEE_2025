import argparse
import os
from collections import defaultdict
import gzip

# récupération des arguments
parser = argparse.ArgumentParser()
parser.add_argument("clean2_Set")
parser.add_argument("sample_evo_scores")
args = parser.parse_args()

bed_file = args.clean2_Set
output_file = args.sample_evo_scores
dir_path = os.path.dirname(output_file)
score_type = os.path.basename(dir_path)
path = "/Users/children/Desktop/Stage_2025/zebrafish/Data"
suffix = f".{score_type}12way.bedGraph.gz"


# Tri du fichier BED -> coordonnées dans ordre croissant
def sorted_dictionary(file):
    dic = defaultdict(list)
    open_func = gzip.open if file.endswith('.gz') else open
    with open_func(file, 'rt') as f:
        for line in f.readlines():
            line = line.strip("\n").split("\t")
            chrom = line[0]
            if file == bed_file:
                pos = (int(line[1]), int(line[2]), str(line[3]))
            else:
                pos = (int(line[1]), int(line[2]), str(line[3]))
            dic[chrom].append(pos)
    if file == bed_file:
        for k in dic.keys():
            dic[k] = list(set(tuple(x) for x in dic[k]))
            dic[k].sort(key=lambda x: x[0])
    return dic

sorted_bed_file = sorted_dictionary(bed_file)

# Overlap interest to reference
dic_output = defaultdict(list)
count_ID_high = 0
count_score_high = 0
for chrom in sorted_bed_file.keys():
    score = f"{path}/danRer_{score_type}_all_chr/{chrom}.danRer11{suffix}"
    if not os.path.exists(score):
        print(f"{chrom} doesn't exist {score}!")
        continue
    else:
        print(chrom)

    score_dic = sorted_dictionary(score)
    first_i = 0
    for pos in sorted_bed_file[chrom]:
        start, end, ID = pos[0], pos[1], str(pos[2])
        i = first_i
        while i < len(score_dic[chrom]) and score_dic[chrom][i][1] <= start:
            i += 1
        first_i = i
        while i < len(score_dic[chrom]) and score_dic[chrom][i][0] < end:
            dic_output[ID].append(score_dic[chrom][i][2])
            i += 1
        ID_len = pos[1]-pos[0]
        score_len = len(dic_output[ID])
        if ID_len > score_len:
            count_ID_high += 1
        if ID_len < score_len:
            count_score_high += 1

with open(output_file, 'w') as f:
    for ID in dic_output.keys():
        f.write(ID + "\t" + '\t'.join(dic_output[ID]) + "\n")
