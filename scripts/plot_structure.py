import sys
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
import numpy as np

# Аргументы от Snakemake
eigenvec_f = sys.argv[1]
eigenval_f = sys.argv[2]
mibs_f = sys.argv[3]
mibs_id_f = sys.argv[4]
out_pca = sys.argv[5]
out_tree = sys.argv[6]

# --- PCA ---
try:
    pca_data = pd.read_csv(eigenvec_f, sep='\s+', header=None)
    eigenval = pd.read_csv(eigenval_f, header=None)
    pca_data.columns = ['FID', 'IID'] + [f'PC{i+1}' for i in range(len(pca_data.columns)-2)]

    total_variance = eigenval[0].sum()
    pc1_var = round((eigenval[0][0] / total_variance) * 100, 2)
    pc2_var = round((eigenval[0][1] / total_variance) * 100, 2)

    plt.figure(figsize=(10, 8))
    plt.scatter(pca_data['PC1'], pca_data['PC2'], s=50)
    for i, txt in enumerate(pca_data['IID']):
        plt.annotate(txt, (pca_data['PC1'][i], pca_data['PC2'][i]), fontsize=8, alpha=0.75)
    plt.xlabel(f'PC1 ({pc1_var}%)'); plt.ylabel(f'PC2 ({pc2_var}%)')
    plt.grid(True)
    plt.savefig(out_pca)
except Exception as e:
    print(f"Error plotting PCA: {e}")

# --- Dendrogram ---
try:
    dist_matrix = np.loadtxt(mibs_f)
    ids_ibs = pd.read_csv(mibs_id_f, sep='\s+', header=None)
    num_samples = len(ids_ibs)
    dissimilarities = dist_matrix[np.tril_indices(num_samples, k=-1)]
    linked = linkage(1 - dissimilarities, method='ward')

    plt.figure(figsize=(12, 7))
    dendrogram(linked, orientation='top', labels=ids_ibs[1].tolist(), distance_sort='descending', show_leaf_counts=True)
    plt.xticks(rotation=90); plt.tight_layout()
    plt.savefig(out_tree)
except Exception as e:
    print(f"Error plotting Tree: {e}")
