import sys
from cyvcf2 import VCF

vcf_path = sys.argv[1]
vcf = VCF(vcf_path)

# Допустим, первые 10 образцов - современные, остальные - исторические
counts = {'modern': 0, 'historical': 0}

for var in vcf:
    # Берем только вредные (LoF/Moderate)
    impact = var.INFO.get('ANN').split('|')[2]
    if impact not in ['HIGH', 'MODERATE']: continue
    
    gts = var.gt_types # 0=HomRef, 1=Het, 2=HomAlt
    # Если мы поляризовали VCF, то Alt = Derived
    counts['modern'] += sum([1 for gt in gts[:10] if gt == 1]) + sum([2 for gt in gts[:10] if gt == 2])
    counts['historical'] += sum([1 for gt in gts[10:] if gt == 1]) + sum([2 for gt in gts[10:] if gt == 2])

rxy = counts['modern'] / counts['historical']
print(f"Rxy (Modern/Historical): {rxy}")
