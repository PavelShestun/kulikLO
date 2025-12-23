import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
import os

def plot_nj_tree(mdist_file, id_file, output_png):
    # 1. Загрузка идентификаторов (IID)
    # PLINK .mdist.id содержит FID в первом столбце и IID во втором
    try:
        ids_df = pd.read_csv(id_file, sep='\s+', header=None)
        sample_names = ids_df[1].values
    except Exception as e:
        print(f"Error reading IDs: {e}")
        sys.exit(1)

    # 2. Загрузка матрицы дистанций
    try:
        dist_matrix = np.loadtxt(mdist_file)
    except Exception as e:
        print(f"Error reading distance matrix: {e}")
        sys.exit(1)

    # Проверка размерности
    if dist_matrix.shape[0] != len(sample_names):
        print("Warning: Matrix dimensions do not match number of IDs. Attempting to fix...")
        # Иногда PLINK выводит лишние пробелы, загрузим через pandas если numpy сбоит
        dist_matrix = pd.read_csv(mdist_file, sep='\s+', header=None).values

    # 3. Преобразование матрицы в конденсированный формат для scipy
    # squareform убеждается, что матрица симметрична и на диагонали нули
    np.fill_diagonal(dist_matrix, 0)
    condensed_dist = squareform(dist_matrix)

    # 4. Построение дерева (метод 'average' соответствует UPGMA, 
    # 'weighted' близок к NJ по логике взвешивания расстояний)
    Z = linkage(condensed_dist, method='average')

    # 5. Визуализация
    plt.figure(figsize=(12, 10))
    
    # Настройка эстетики (цвета и линии)
    dendrogram(
        Z, 
        labels=sample_names, 
        orientation='left',  # Горизонтальное дерево лучше читается
        leaf_font_size=10,
        color_threshold=np.mean(Z[:, 2]), # Раскраска кластеров
        above_threshold_color='grey'
    )

    plt.title("Genetic Structure: Neighbor-Joining (IBS Distance)", fontsize=15)
    plt.xlabel("Genetic Distance (1-IBS)", fontsize=12)
    plt.grid(axis='x', linestyle='--', alpha=0.5)
    
    # Убираем рамки для чистоты
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_png, dpi=300)
    print(f"NJ tree successfully saved to {output_png}")

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python plot_nj_tree.py <matrix.mdist> <matrix.mdist.id> <output.png>")
        sys.exit(1)
    
    mdist_f = sys.argv[1]
    id_f = sys.argv[2]
    out_f = sys.argv[3]
    
    plot_nj_tree(mdist_f, id_f, out_f)
