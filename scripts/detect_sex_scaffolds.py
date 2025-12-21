import sys
import pandas as pd

# Аргументы: 1=blast_Z, 2=blast_W, 3=output_list
blast_z_file = sys.argv[1]
blast_w_file = sys.argv[2]
output_file = sys.argv[3]

def get_top_scaffolds(blast_file, top_n=5):
    try:
        # Читаем формат outfmt 6
        df = pd.read_csv(blast_file, sep='\t', header=None, 
                         names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
                                'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
        # Суммируем длину совпадений для каждого скаффолда
        scaffold_counts = df.groupby('sseqid')['length'].sum().sort_values(ascending=False)
        return scaffold_counts.head(top_n).index.tolist()
    except pd.errors.EmptyDataError:
        return []

z_scaffolds = get_top_scaffolds(blast_z_file)
w_scaffolds = get_top_scaffolds(blast_w_file)

# Объединяем уникальные
sex_scaffolds = sorted(list(set(z_scaffolds + w_scaffolds)))

with open(output_file, 'w') as f:
    for scaf in sex_scaffolds:
        f.write(f"{scaf}\n")

print(f"Detected {len(sex_scaffolds)} sex scaffolds: {sex_scaffolds}")
