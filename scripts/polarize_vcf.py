import sys
from cyvcf2 import VCF, Writer
import numpy as np

def polarize(vcf_in, outgroup_vcf, vcf_out):
    # Открываем основной VCF и VCF внешней группы
    # Внешняя группа ДОЛЖНА быть индексирована (tabix)
    target_vcf = VCF(vcf_in)
    out_vcf = VCF(outgroup_vcf)
    
    # Добавляем информацию в хедер о поляризации
    target_vcf.add_info_to_header({
        'ID': 'AA', 'Number': '1', 'Type': 'String', 
        'Description': 'Ancestral Allele determined by outgroup'
    })
    
    # Настраиваем запись
    w = Writer(vcf_out, target_vcf)
    
    count_swapped = 0
    count_total = 0
    count_dropped = 0

    print(f"Starting polarization of {vcf_in} using {outgroup_vcf}...")

    for variant in target_vcf:
        count_total += 1
        
        # Ищем эту же позицию в VCF внешней группы
        region = f"{variant.CHROM}:{variant.POS}-{variant.POS}"
        ancestral_allele = None
        
        # Получаем записи из внешней группы для этой позиции
        for out_var in out_vcf(region):
            # Проверяем только если позиции совпадают точно
            if out_var.POS == variant.POS:
                # Берем наиболее частый аллель во внешней группе
                # (в идеале внешняя группа должна быть гомозиготна)
                # gt_types: 0=HOM_REF, 1=HET, 2=HOM_ALT, 3=UNKNOWN
                gts = out_var.gt_types
                if np.mean(gts == 0) > 0.8: # Почти все HOM_REF
                    ancestral_allele = out_var.REF
                elif np.mean(gts == 2) > 0.8: # Почти все HOM_ALT
                    ancestral_allele = out_var.ALT[0]
                break
        
        # Если не нашли аллель во внешней группе - пропускаем сайт (или оставляем как есть)
        # В статье по какапо такие сайты часто исключали для чистоты Rxy
        if ancestral_allele is None:
            count_dropped += 1
            continue

        # ЛОГИКА ПОЛЯРИЗАЦИИ:
        # Если Ancestral == нашему ALT, значит наш REF на самом деле является мутацией (Derived).
        # Нужно поменять их местами.
        
        if ancestral_allele == variant.ALT[0]:
            # 1. Меняем REF и ALT местами
            old_ref = variant.REF
            old_alt = variant.ALT[0]
            variant.REF = old_alt
            variant.ALT = [old_ref]
            
            # 2. Инвертируем генотипы у всех особей
            # 0 -> 1 (был REF, стал ALT)
            # 1 -> 1 (гетерозигота остается гетерозиготой)
            # 2 -> 0 (был ALT, стал REF)
            # Используем прямое манипулирование массивом генотипов cyvcf2
            new_gts = variant.genotypes
            for i in range(len(new_gts)):
                gt = new_gts[i] # Формат [allele1, allele2, is_phased]
                for j in range(2):
                    if gt[j] == 0: gt[j] = 1
                    elif gt[j] == 1: gt[j] = 0
                new_gts[i] = gt
            variant.genotypes = new_gts
            
            count_swapped += 1
        
        # Если Ancestral == нашему REF, ничего менять не надо, ALT и так Derived.
        # Если Ancestral не совпал ни с чем (третий аллель), сайт обычно выбрасывают.
        elif ancestral_allele != variant.REF:
            count_dropped += 1
            continue
            
        # Добавляем метку предкового аллеля в INFO
        variant.INFO["AA"] = ancestral_allele
        w.write_record(variant)

    w.close()
    target_vcf.close()
    print(f"Done!")
    print(f"Total processed: {count_total}")
    print(f"Swapped (REF/ALT flip): {count_swapped}")
    print(f"Dropped (no outgroup data or mismatched): {count_dropped}")

if __name__ == "__main__":
    vcf_in = sys.argv[1]
    outgroup_vcf = sys.argv[2]
    vcf_out = sys.argv[3]
    polarize(vcf_in, outgroup_vcf, vcf_out)
