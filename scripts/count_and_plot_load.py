import vcf
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import numpy as np

VCF_FILE = sys.argv[1]
OUT_PLOT = sys.argv[2]

try:
    vcf_reader = vcf.Reader(filename=VCF_FILE)
except:
    print("Error opening VCF or empty file.")
    sys.exit(0)

load_counts = {sample: {'syn_het': 0, 'syn_hom': 0, 'mis_het': 0, 'mis_hom': 0, 'lof_het': 0, 'lof_hom': 0} for sample in vcf_reader.samples}
count = 0

for record in vcf_reader:
    count += 1
    if 'ANN' not in record.INFO: continue
    ann_str = record.INFO['ANN'][0]
    try: impact = ann_str.split('|')[2]
    except IndexError: continue
    
    category = None
    if impact == 'LOW': category = 'syn'
    elif impact == 'MODERATE': category = 'mis'
    elif impact == 'HIGH': category = 'lof'
    if category is None: continue

    for sample in record.samples:
        if not sample.called: continue
        if sample.gt_type == 1: load_counts[sample.sample][f'{category}_het'] += 1
        elif sample.gt_type == 2: load_counts[sample.sample][f'{category}_hom'] += 1

if count == 0:
    print("No variants found in annotated VCF.")
    plt.figure()
    plt.text(0.5, 0.5, "No Data found")
    plt.savefig(OUT_PLOT)
    sys.exit(0)

df = pd.DataFrame.from_dict(load_counts, orient='index')
df.index.name = 'IID'; df.reset_index(inplace=True)
df['syn_total'] = df['syn_hom'] + df['syn_het']
df['mis_total'] = df['mis_hom'] + df['mis_het']
df['lof_total'] = df['lof_hom'] + df['lof_het']

# Защита от деления на ноль
df['Ratio_LoF_Syn'] = df.apply(lambda row: row['lof_total'] / row['syn_total'] if row['syn_total'] > 0 else 0, axis=1)
df['Ratio_Mis_Syn'] = df.apply(lambda row: row['mis_total'] / row['syn_total'] if row['syn_total'] > 0 else 0, axis=1)

plot_df = df[['Ratio_LoF_Syn', 'Ratio_Mis_Syn']].melt(var_name='Type', value_name='Ratio')

plt.figure(figsize=(10, 8))
if not plot_df.empty:
    # Исправление для Seaborn: hue=Type и legend=False
    sns.boxplot(x='Type', y='Ratio', hue='Type', data=plot_df, palette="Pastel1", legend=False)
    sns.stripplot(x='Type', y='Ratio', data=plot_df, color='black', alpha=0.6, jitter=True)

plt.savefig(OUT_PLOT)
