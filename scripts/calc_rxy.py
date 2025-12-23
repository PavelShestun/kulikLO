import sys
from cyvcf2 import VCF

vcf_path = sys.argv[1]
out_path = sys.argv[2]
vcf = VCF(vcf_path)
samples = vcf.samples

# Словарь для хранения результатов по каждой особи
# Считаем количество Derived (мутантных) аллелей
results = {s: {'lof': 0, 'missense': 0} for s in samples}

print(f"Calculating mutational load for {len(samples)} samples...")

for var in vcf:
    ann = var.INFO.get('ANN')
    if not ann: continue
    
    # Извлекаем тип воздействия мутации
    impact = ann.split('|')[2]
    
    # Логика: считаем только вредные мутации (LoF и Moderate)
    category = None
    if impact == 'HIGH': category = 'lof'
    elif impact == 'MODERATE': category = 'missense'
    
    if category:
        gts = var.gt_types # 0=HomRef, 1=Het, 2=HomAlt
        for i, gt in enumerate(gts):
            if gt == 1: # Гетерозигота (1 копия мутации)
                results[samples[i]][category] += 1
            elif gt == 2: # Гомозигота (2 копии мутации)
                results[samples[i]][category] += 2

# Сохраняем результат в таблицу
with open(out_path, 'w') as f:
    f.write("SampleID\tLoF_count\tMissense_count\n")
    for s in samples:
        f.write(f"{s}\t{results[s]['lof']}\t{results[s]['missense']}\n")

print(f"Done! Statistics saved to {out_path}")
