import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

roh_file = sys.argv[1]
out_plot = sys.argv[2]

# Проверка, пустой ли файл
if os.path.getsize(roh_file) == 0:
    print("ROH file is empty. No ROH detected.")
    plt.figure()
    plt.text(0.5, 0.5, "No ROH Detected (Data Empty)", ha='center')
    plt.savefig(out_plot)
    sys.exit()

try:
    df = pd.read_csv(roh_file, sep='\s+')
    if df.empty:
        raise ValueError("Empty DataFrame")
        
    # ОЧИСТКА ИМЕН: убираем путь и расширение, оставляем только ID
    df['IID'] = df['IID'].astype(str).apply(lambda x: os.path.basename(x).split('.')[0])
except Exception as e:
    print(f"Error parsing ROH file: {e}")
    plt.figure()
    plt.text(0.5, 0.5, f"No ROH segments found\nCheck filtering parameters", ha='center')
    plt.savefig(out_plot)
    sys.exit()

GENOME_LEN = 1.2e9 
# Считаем сумму ROH (KB) для каждой особи
froh = df.groupby('IID')['KB'].sum() * 1000 / GENOME_LEN
froh = froh.reset_index(name='F_ROH')
n_roh = df.groupby('IID')['KB'].count().reset_index(name='N_ROH')

data = pd.merge(froh, n_roh, on='IID')

plt.figure(figsize=(12, 7))
if not data.empty:
    sns.scatterplot(data=data, x='N_ROH', y='F_ROH', s=200, color='darkred', alpha=0.7, edgecolor='black')

    for i in range(data.shape[0]):
        plt.text(data.N_ROH[i] + 0.1, data.F_ROH[i], data.IID[i], fontsize=9, verticalalignment='bottom')

plt.xlabel("Number of ROH segments (N_ROH)")
plt.ylabel("Inbreeding Coefficient (F_ROH)")
plt.title("Genomic Inbreeding: Runs of Homozygosity Analysis")
plt.grid(True, alpha=0.3, linestyle='--')
plt.savefig(out_plot, dpi=300)
