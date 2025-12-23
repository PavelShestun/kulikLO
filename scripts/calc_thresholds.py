import sys
import numpy as np
from cyvcf2 import VCF

# Использование: python scripts/calc_thresholds.py resources/merged.vcf
vcf_file = sys.argv[1]

print(f"Reading {vcf_file} to calculate depth stats...")
depths = []
vcf = VCF(vcf_file)

# Читаем первые 100,000 вариантов для скорости
for i, variant in enumerate(vcf):
    if i > 100000: break
    dp = variant.INFO.get('DP')
    if dp:
        depths.append(dp)

depths = np.array(depths)
mean_dp = np.mean(depths)
std_dp = np.std(depths)

print("\n" + "="*30)
print(f"  MEAN DEPTH (Global): {mean_dp:.2f}")
print("="*30)
print("Recommended config settings:")
print(f"  min_depth: {int(mean_dp * 0.3)}  (Mean * 0.3)")
print(f"  max_depth: {int(mean_dp * 1.8)}  (Mean * 1.8 - cutoff for paralogs)")
print("="*30 + "\n")
