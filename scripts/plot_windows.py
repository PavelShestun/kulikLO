import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

pi_file = sys.argv[1]
tajima_file = sys.argv[2]
out_plot = sys.argv[3]

# Читаем Pi
try:
    df_pi = pd.read_csv(pi_file, sep='\t')
    df_taj = pd.read_csv(tajima_file, sep='\t')
except:
    plt.figure(); plt.savefig(out_plot); sys.exit()

fig, axes = plt.subplots(2, 1, figsize=(12, 10))

# График 1: Pi
# Фильтруем слишком большие значения (ошибки) и NaN
df_pi = df_pi[df_pi['PI'] < 0.1].dropna()
sns.histplot(df_pi['PI'], bins=50, color='teal', ax=axes[0], kde=True)
axes[0].set_title('Nucleotide Diversity (Pi) Distribution')
axes[0].set_xlabel('Pi (per site)')

# График 2: Tajima's D
df_taj = df_taj.dropna()
sns.histplot(df_taj['TajimaD'], bins=50, color='purple', ax=axes[1], kde=True)
axes[1].set_title("Tajima's D Distribution")
axes[1].set_xlabel("Tajima's D")
axes[1].axvline(x=0, color='black', linestyle='--')

plt.tight_layout()
plt.savefig(out_plot)
