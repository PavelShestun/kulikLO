import subprocess
import os

configfile: "config/config.yaml"

# Пытаемся получить список образцов напрямую из VCF файла
try:
    VCF_PATH = config["resources"]["vcf"]
    # Если файл существует, читаем образцы. Если нет - пустой список (Snakemake сам разберется при запуске)
    if os.path.exists(VCF_PATH):
        SAMPLES = subprocess.check_output(f"bcftools query -l {VCF_PATH}", shell=True).decode().split()
    else:
        SAMPLES = []
except:
    SAMPLES = []

rule all:
    input:
        "results/qc/clean_strict.vcf.gz.tbi",    # Подтверждение окончания фильтрации
        "results/plots/PCA_plot.png",            # График PCA
        "results/plots/NJ_tree.png",             # Дерево
        "results/plots/sNMF_entropy.png",        # Структура (sNMF)
        "results/plots/ROH_Inbreeding_Final.png",# Инбридинг
        "results/load/Rxy_stats.txt",            # Мутационный груз
        "results/plots/Genomic_Windows.png",     # Pi и Tajima's D
        "results/dadi/best_model_plot.png",      # Демография dadi
        "results/models/slim_purging_results.png"# Моделирование SLiM

# --- 1. Определение пола (Blast по Z/W) ---
rule download_sex_refs:
    output: z="resources/ref_sex/chrZ.fna", w="resources/ref_sex/chrW.fna"
    params: z_url=config["resources"]["ref_Z_url"], w_url=config["resources"]["ref_W_url"]
    shell: 
        "wget -O {output.z}.gz {params.z_url} && gunzip -f {output.z}.gz; "
        "wget -O {output.w}.gz {params.w_url} && gunzip -f {output.w}.gz"

rule blast_db:
    input: config["resources"]["genome"]
    output: multiext("resources/blast_db/calidris", ".ndb", ".nhr", ".nsq")
    shell: "makeblastdb -in {input} -dbtype nucl -out resources/blast_db/calidris"

rule blast_search:
    input: z="resources/ref_sex/chrZ.fna", w="resources/ref_sex/chrW.fna", db=multiext("resources/blast_db/calidris", ".ndb", ".nhr", ".nsq")
    output: z_res="results/sex_check/blast_Z.txt", w_res="results/sex_check/blast_W.txt"
    threads: 4
    shell: 
        "blastn -query {input.z} -db resources/blast_db/calidris -out {output.z_res} -outfmt 6 -evalue 1e-10 -num_threads {threads}; "
        "blastn -query {input.w} -db resources/blast_db/calidris -out {output.w_res} -outfmt 6 -evalue 1e-10 -num_threads {threads}"

rule identify_sex_scaffolds:
    input: z="results/sex_check/blast_Z.txt", w="results/sex_check/blast_W.txt"
    output: "results/sex_check/sex_scaffolds.txt"
    shell: "python scripts/detect_sex_scaffolds.py {input.z} {input.w} {output}"

# --- 2. ЖЕСТКАЯ ФИЛЬТРАЦИЯ ---
rule filter_vcf_strict:
    input: 
        vcf=config["resources"]["vcf"],
        sex="results/sex_check/sex_scaffolds.txt"
    output: 
        vcf="results/qc/clean_strict.vcf.gz",
        tbi="results/qc/clean_strict.vcf.gz.tbi"
    params:
        min_dp=config["filtering"]["min_depth"],
        max_dp=config["filtering"]["max_depth"],
        ab_min=config["filtering"]["min_allele_balance"],
        snp_gap=config["filtering"]["snp_gap"],
        hwe=config["filtering"]["hwe_threshold"]
    shell:
        """
        # Убрали python скрипт из цепочки, чтобы избежать ошибок потока
        bcftools view -m2 -M2 -v snps {input.vcf} | \
        bcftools filter -e 'QUAL < 30 || INFO/DP < {params.min_dp} || INFO/DP > {params.max_dp} || HWE < {params.hwe}' | \
        bcftools filter -g {params.snp_gap} | \
        bcftools view -O z -o {output.vcf}
        
        bcftools index -t {output.vcf}
        """

# --- 3. Структура популяции (PCA + NJ Tree + sNMF) ---
rule ld_pruning:
    input: "results/qc/clean_strict.vcf.gz"
    output: 
        vcf="results/qc/pruned.vcf.gz",
        prune_in="results/qc/pruning.prune.in"
    shell:
        """
        # --double-id лечит ошибку с именами образцов
        # --vcf-half-call лечит ошибки некорректных генотипов
        plink --vcf {input} --allow-extra-chr --double-id --vcf-half-call missing \
              --indep-pairwise 50 10 0.2 --out results/qc/pruning
              
        plink --vcf {input} --allow-extra-chr --double-id --vcf-half-call missing \
              --extract results/qc/pruning.prune.in --recode vcf --out results/qc/pruned
              
        bgzip -f results/qc/pruned.vcf && bcftools index -t results/qc/pruned.vcf.gz
        """

rule pca_analysis:
    input: "results/qc/pruned.vcf.gz"
    output: vec="results/qc/pca.eigenvec", val="results/qc/pca.eigenval"
    shell: "plink --vcf {input} --allow-extra-chr --double-id --pca --out results/qc/pca"

# Рассчитываем матрицу дистанций для дерева
rule calc_dist_matrix:
    input: "results/qc/pruned.vcf.gz"
    output: mdist="results/qc/dist_matrix.mdist", id="results/qc/dist_matrix.mdist.id"
    shell: "plink --vcf {input} --allow-extra-chr --double-id --distance 1-ibs square --out results/qc/dist_matrix"

rule nj_tree:
    input: mdist="results/qc/dist_matrix.mdist", id="results/qc/dist_matrix.mdist.id"
    output: "results/plots/NJ_tree.png"
    shell: "python scripts/plot_nj_tree.py {input.mdist} {input.id} {output}"

# Добавленное правило: Рисуем PCA (использует скрипт plot_structure.py)
rule plot_pca_structure:
    input: 
        vec="results/qc/pca.eigenvec", 
        val="results/qc/pca.eigenval",
        mdist="results/qc/dist_matrix.mdist", 
        id="results/qc/dist_matrix.mdist.id"
    output: 
        pca_plot="results/plots/PCA_plot.png",
        tree_plot_dummy="results/plots/Tree_structure_dummy.png" # Скрипт делает два графика сразу
    shell:
        "python scripts/plot_structure.py {input.vec} {input.val} {input.mdist} {input.id} {output.pca_plot} {output.tree_plot_dummy}"

rule run_snmf:
    input: "results/qc/pruned.vcf.gz"
    output: "results/plots/sNMF_entropy.png"
    params: 
        k_min=config["snmf"]["K_min"], 
        k_max=config["snmf"]["K_max"], 
        reps=config["snmf"]["repetitions"]
    shell:
        """
        # Очистка старых проектов sNMF (они могут мешать)
        rm -rf results/qc/pruned_for_lea.vcf.snmf/
        rm -f results/qc/pruned_for_lea.vcf.geno
        
        zcat {input} > results/qc/pruned_for_lea.vcf
        Rscript scripts/run_snmf.R results/qc/pruned_for_lea.vcf {params.k_min} {params.k_max} {params.reps} {output}
        rm results/qc/pruned_for_lea.vcf
        """

# --- 4. Демографическая история (dadi) ---
rule dadi_inference:
    input: "results/qc/pruned.vcf.gz"
    output: params="results/dadi/best_params.txt", plot="results/dadi/best_model_plot.png"
    shell: 
        """
        # Генерируем правильный pop_map: берем список образцов из VCF 
        # и добавляем к каждому колонку "Spoonbill"
        bcftools query -l {input} | awk '{{print $1"\tSpoonbill"}}' > pop_map_real.txt
        
        # Запускаем dadi
        python scripts/dadi_inference.py {input} pop_map_real.txt {output.params} {output.plot}
        
        # Удаляем временный файл
        rm pop_map_real.txt
        """

# --- 5. Мутационный груз (SnpEff + Rxy) ---
rule build_snpeff_db:
    input: g=config["resources"]["genome"], a=config["resources"]["genes"]
    output: directory("results/snpeff_data/calidris")
    shell:
        """
        mkdir -p results/snpeff_data/calidris
        cp {input.g} results/snpeff_data/calidris/sequences.fa
        cp {input.a} results/snpeff_data/calidris/genes.gff
        
        # --- ИСПРАВЛЕНИЕ ---
        # 1. Указываем SnpEff, что папка с данными - это results/snpeff_data
        echo "data.dir = results/snpeff_data/" > snpEff_local.config
        # 2. Добавляем название генома
        echo "calidris.genome : Calidris" >> snpEff_local.config
        
        # Строим базу
        snpEff build -c snpEff_local.config -gff3 -v calidris -noCheckCds -noCheckProtein
        """

# --- ПОЛЯРИЗАЦИЯ  ---
# Создаем VCF для внешней группы путем выравнивания её генома на референс
rule map_outgroup:
    input:
        ref=config["resources"]["genome"],
        out_gen=config["resources"]["outgroup_genome"]
    output:
        bam="results/load/outgroup_mapped.bam"
    threads: 1
    shell:
        """
        # 1. Выравнивание в несжатый BAM (minimap2 съест ~7ГБ и закроется)
        minimap2 -K 20M -t {threads} -ax asm5 {input.ref} {input.out_gen} | \
        samtools view -b - > results/load/unsorted.bam
        
        # 2. Сортировка (теперь памяти будет много, выделяем 4ГБ)
        samtools sort -m 4G -@ {threads} -o {output.bam} results/load/unsorted.bam
        
        # 3. Удаление временного файла
        rm results/load/unsorted.bam
        """

rule call_outgroup_vcf:
    input:
        ref=config["resources"]["genome"],
        bam="results/load/outgroup_mapped.bam"
    output:
        vcf=config["resources"]["outgroup_vcf"],
        tbi=config["resources"]["outgroup_vcf"] + ".tbi"
    shell:
        """
        bcftools mpileup -f {input.ref} {input.bam} | \
        bcftools call -mv -Oz -o {output.vcf}
        bcftools index -t {output.vcf}
        """

rule polarize_vcf:
    input: 
        vcf="results/qc/clean_strict.vcf.gz",
        outgroup=config["resources"]["outgroup_vcf"] # Теперь это созданный выше файл
    output: 
        vcf="results/load/polarized.vcf.gz",
        tbi="results/load/polarized.vcf.gz.tbi"
    shell: 
        """
        python scripts/polarize_vcf.py {input.vcf} {input.outgroup} {output.vcf}
        bcftools index -t {output.vcf}
        """

rule annotate_load:
    input: 
        vcf="results/load/polarized.vcf.gz", # Возвращаем polarized
        db="results/snpeff_data/calidris"
    output: "results/load/annotated.vcf.gz"
    shell: "snpEff ann -c snpEff_local.config -v calidris {input.vcf} | bgzip > {output} && bcftools index -t {output}"

rule calculate_rxy:
    input: "results/load/annotated.vcf.gz"
    output: "results/load/Rxy_stats.txt"
    shell: "python scripts/calc_rxy.py {input} {output}"

# --- 6. Инбридинг (ROH) ---
rule calc_roh_final:
    input: "results/qc/clean_strict.vcf.gz"
    output: "results/plots/ROH_Inbreeding_Final.png"
    params:
        kb=config["inbreeding_roh"]["min_kb"],
        snps=config["inbreeding_roh"]["min_snps"]
    shell:
        """
        plink --vcf {input} --allow-extra-chr --double-id --vcf-half-call missing \
              --homozyg --homozyg-kb {params.kb} --homozyg-snp {params.snps} \
              --homozyg-window-het 1 --homozyg-density 50 \
              --out results/qc/roh_final
        
        if [ -f results/qc/roh_final.hom ]; then
            python scripts/plot_roh.py results/qc/roh_final.hom {output}
        else
            # Если ROH не найдены (редко, но бывает), создаем пустой файл
            touch {output}
        fi
        """

# --- 7. Геномное сканирование ---
rule genomic_windows:
    input: "results/qc/clean_strict.vcf.gz"
    output: "results/plots/Genomic_Windows.png"
    params: w=config["window_stats"]["window_size"]
    shell:
        """
        vcftools --gzvcf {input} --window-pi {params.w} --out results/qc/windows
        vcftools --gzvcf {input} --TajimaD {params.w} --out results/qc/windows
        python scripts/plot_windows.py results/qc/windows.windowed.pi results/qc/windows.Tajima.D {output}
        """

# --- 8. SLiM моделирование (Purging) ---
rule slim_purging:
    input: "results/dadi/best_params.txt"
    output: "results/models/slim_purging_results.png"
    # Используем run_final_slim.py (или run_slim_dynamic.py, в зависимости от того, какой скрипт вы создали)
    shell: "python scripts/run_final_slim.py {output}"
