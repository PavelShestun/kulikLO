import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from cyvcf2 import VCF
import numpy as np

VCF_FILE = sys.argv[1]
OUT_PLOT = sys.argv[2]

try:
    vcf = VCF(VCF_FILE)
    samples = vcf.samples
except Exception as e:
    print(f"Error opening VCF: {e}")
    sys.exit(1)

# Инициализация счетчиков
# load_counts = {'SampleName': {'syn': 0, 'mis': 0, 'lof': 0}}
load_counts = {s: {'syn': 0, 'mis': 0, 'lof': 0} for s in samples}

print("Parsing VCF...")
for variant in vcf:
    # Пропускаем, если нет аннотации
    ann_info = variant.INFO.get('ANN')
    if not ann_info:
        continue
    
    # Берем первую аннотацию (обычно канонический транскрипт)
    # Формат SnpEff: Allele|Annotation|Impact|GeneName...
    ann_parts = ann_info.split(',')[0].split('|')
    if len(ann_parts) < 3:
        continue
        
    impact = ann_parts[2]
    
    category = None
    if impact == 'LOW': category = 'syn'
    elif impact == 'MODERATE': category = 'mis'
    elif impact == 'HIGH': category = 'lof'
    
    if category is None:
        continue

    # Итерация по генотипам (намного быстрее в cyvcf2)
    # variant.gt_types: 0=HOM_REF, 1=HET, 2=HOM_ALT, 3=UNKNOWN
    for i, gt in enumerate(variant.gt_types):
        if gt == 1 or gt == 2: # Если есть альтернативный аллель (HET или HOM)
            # В данном простом анализе считаем "нагрузкой" наличие аллеля
            # Можно усложнить: HOM = 2, HET = 1
            load_counts[samples[i]][category] += 1

# Создание DataFrame
df = pd.DataFrame.from_dict(load_counts, orient='index')
df.index.name = 'IID'
df.reset_index(inplace=True)

# Расчет отношений
# Добавляем малый эпсилон, чтобы не делить на ноль
df['Ratio_LoF_Syn'] = df['lof'] / (df['syn'] + 1e-9)
df['Ratio_Mis_Syn'] = df['mis'] / (df['syn'] + 1e-9)

plot_df = df[['Ratio_LoF_Syn', 'Ratio_Mis_Syn']].melt(var_name='Type', value_name='Ratio')

plt.figure(figsize=(10, 8))
if not plot_df.empty:
    sns.boxplot(x='Type', y='Ratio', hue='Type', data=plot_df, palette="Pastel1", legend=False)
    sns.stripplot(x='Type', y='Ratio', data=plot_df, color='black', alpha=0.6, jitter=True)
    plt.title("Mutational Load Ratios (Relative to Synonymous)")
else:
    plt.text(0.5, 0.5, "No variants found")

plt.savefig(OUT_PLOT)
