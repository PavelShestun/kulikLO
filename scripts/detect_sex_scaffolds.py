import sys
import pandas as pd

blast_z_file = sys.argv[1]
blast_w_file = sys.argv[2]
output_file = sys.argv[3]

def get_top_scaffolds(blast_file, top_n=10):
    try:
        # Проверяем, пустой ли файл
        import os
        if os.path.getsize(blast_file) == 0:
            return []
            
        df = pd.read_csv(blast_file, sep='\t', header=None, 
                         names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
                                'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
        
        # Фильтруем слабые совпадения
        df = df[(df['pident'] > 80) & (df['length'] > 500)]
        
        if df.empty:
            return []

        # Суммируем длину выравниваний
        scaffold_counts = df.groupby('sseqid')['length'].sum().sort_values(ascending=False)
        return scaffold_counts.head(top_n).index.tolist()
    except Exception as e:
        print(f"Warning processing {blast_file}: {e}")
        return []

z_scaffolds = get_top_scaffolds(blast_z_file)
w_scaffolds = get_top_scaffolds(blast_w_file)

sex_scaffolds = sorted(list(set(z_scaffolds + w_scaffolds)))

with open(output_file, 'w') as f:
    for scaf in sex_scaffolds:
        f.write(f"{scaf}\n")

print(f"Detected {len(sex_scaffolds)} sex scaffolds.")
