configfile: "config/config.yaml"

rule all:
    input:
        "results/qc/stats_final.txt",
        "results/plots/PCA_plot.png",
        "results/plots/sNMF_entropy.png",
        "results/dadi/best_model_plot.png",
        "results/load/Mutational_Load_plot_final.png",
        "results/models/slim_dynamic_result.png",
        "results/plots/ROH_Inbreeding.png",
        "results/plots/Genomic_Windows.png"

# --- 1. Определение пола ---
rule download_sex_refs:
    output: z="resources/ref_sex/chrZ.fna", w="resources/ref_sex/chrW.fna"
    params: z_url=config["resources"]["ref_Z_url"], w_url=config["resources"]["ref_W_url"]
    shell: "wget -O {output.z}.gz {params.z_url} && gunzip {output.z}.gz; wget -O {output.w}.gz {params.w_url} && gunzip {output.w}.gz"

rule blast_db:
    input: config["resources"]["genome"]
    output: multiext("resources/blast_db/calidris", ".ndb", ".nhr", ".nsq")
    shell: "makeblastdb -in {input} -dbtype nucl -out resources/blast_db/calidris"

rule blast_search:
    input: z="resources/ref_sex/chrZ.fna", w="resources/ref_sex/chrW.fna", db=multiext("resources/blast_db/calidris", ".ndb", ".nhr", ".nsq")
    output: z_res="results/sex_check/blast_Z.txt", w_res="results/sex_check/blast_W.txt"
    threads: 8
    shell: "blastn -query {input.z} -db resources/blast_db/calidris -out {output.z_res} -outfmt 6 -evalue 1e-10 -num_threads {threads}; blastn -query {input.w} -db resources/blast_db/calidris -out {output.w_res} -outfmt 6 -evalue 1e-10 -num_threads {threads}"

rule identify_sex_scaffolds:
    input: z="results/sex_check/blast_Z.txt", w="results/sex_check/blast_W.txt"
    output: "results/sex_check/sex_scaffolds.txt"
    shell: "python scripts/detect_sex_scaffolds.py {input.z} {input.w} {output}"

# --- 2. Фильтрация ---
rule filter_vcf_full:
    input: vcf=config["resources"]["vcf"], sex="results/sex_check/sex_scaffolds.txt"
    output: "results/qc/clean_full.vcf.gz"
    shell: "bcftools view {input.vcf} -O z -o {output} && bcftools index {output}"

rule qc_stats:
    input: "results/qc/clean_full.vcf.gz"
    output: "results/qc/stats_final.txt"
    shell: "bcftools stats {input} > {output}"

# --- 3. PCA (Исправлено рисование) ---
rule prune_data:
    input: "results/qc/clean_full.vcf.gz"
    output: pruned_vcf="results/qc/pruned.vcf", vec="results/qc/pca.eigenvec", val="results/qc/pca.eigenval"
    shell:
        """
        bcftools annotate --set-id '%CHROM:%POS' {input} -O z -o results/qc/annotated.vcf.gz
        vcftools --gzvcf results/qc/annotated.vcf.gz --plink --out results/qc/temp_plink
        sed -i 's/^0/1/' results/qc/temp_plink.map
        plink --file results/qc/temp_plink --allow-extra-chr --indep-pairwise 50 10 0.2 --out results/qc/pruning
        plink --file results/qc/temp_plink --allow-extra-chr --extract results/qc/pruning.prune.in --make-bed --out results/qc/for_pca
        plink --bfile results/qc/for_pca --recode vcf --out results/qc/pruned
        plink --bfile results/qc/for_pca --allow-extra-chr --pca --out results/qc/pca
        mv results/qc/pca.eigenvec {output.vec}
        mv results/qc/pca.eigenval {output.val}
        rm results/qc/annotated.vcf.gz results/qc/temp_plink*
        """

rule run_pca_plot:
    input: vec="results/qc/pca.eigenvec", val="results/qc/pca.eigenval"
    output: "results/plots/PCA_plot.png"
    shell:
        """
        python -c "
import pandas as pd
import matplotlib.pyplot as plt
import os
pca = pd.read_csv('{input.vec}', sep='\s+', header=None)
pca[1] = pca[1].apply(lambda x: os.path.basename(x).split('.')[0])
plt.figure(figsize=(10,8))
plt.scatter(pca[2], pca[3], s=150, color='blue', alpha=0.6, edgecolor='k')
for i, txt in enumerate(pca[1]):
    plt.annotate(txt, (pca.iloc[i, 2], pca.iloc[i, 3]), fontsize=8)
plt.xlabel('PC1'); plt.ylabel('PC2'); plt.title('PCA Analysis')
plt.grid(True, alpha=0.2)
plt.savefig('{output}')"
        """

# --- 4-8. Остальные анализы ---
rule run_snmf:
    input: "results/qc/pruned.vcf"
    output: plot="results/plots/sNMF_entropy.png"
    params: k_min=config["snmf"]["K_min"], k_max=config["snmf"]["K_max"], reps=config["snmf"]["repetitions"]
    shell: "Rscript scripts/run_snmf.R {input} {params.k_min} {params.k_max} {params.reps} {output.plot} || touch {output.plot}"

rule dadi_inference:
    input: "results/qc/pruned.vcf"
    output: params="results/dadi/best_params.txt", plot="results/dadi/best_model_plot.png"
    shell:
        """
        bcftools query -l {input} > results/dadi/pop.txt
        awk '{{print $1 "\\tSpoonbill"}}' results/dadi/pop.txt > results/dadi/pop_map.txt
        python scripts/dadi_inference.py {input} results/dadi/pop_map.txt {output.params} {output.plot}
        """

rule snpeff_db:
    input: g=config["resources"]["genome"], a=config["resources"]["genes"]
    output: directory("results/snpeff_data/calidris_pygmaea")
    shell: "mkdir -p results/snpeff_data/calidris_pygmaea; cp {input.g} results/snpeff_data/calidris_pygmaea/sequences.fa; cp {input.a} results/snpeff_data/calidris_pygmaea/genes.gff; echo 'data.dir = ./results/snpeff_data/' > snpEff_local.config; echo 'calidris_pygmaea.genome : Calidris pygmaea' >> snpEff_local.config; snpEff build -c snpEff_local.config -gff3 -v -noCheckCds -noCheckProtein calidris_pygmaea"

rule annotate_load:
    input: vcf="results/qc/clean_full.vcf.gz", db="results/snpeff_data/calidris_pygmaea"
    output: "results/load/annotated.vcf"
    shell: "snpEff ann -c snpEff_local.config -v calidris_pygmaea {input.vcf} > {output}"

rule plot_load:
    input: "results/load/annotated.vcf"
    output: "results/load/Mutational_Load_plot_final.png"
    shell: "python scripts/count_and_plot_load.py {input} {output}"

rule slim_dynamic:
    input: "results/dadi/best_params.txt"
    output: "results/models/slim_dynamic_result.png"
    shell: "python scripts/run_slim_dynamic.py {input} {output}"

rule calc_roh:
    input: "results/qc/clean_full.vcf.gz"
    output: hom="results/qc/plink.hom", plot="results/plots/ROH_Inbreeding.png"
    shell:
        """
        bcftools annotate --set-id '%CHROM:%POS' {input} -O z -o results/qc/ann_roh.vcf.gz
        vcftools --gzvcf results/qc/ann_roh.vcf.gz --plink --out results/qc/temp_roh
        sed -i 's/^0/1/' results/qc/temp_roh.map
        plink --file results/qc/temp_roh --allow-extra-chr --homozyg --homozyg-kb 100 --homozyg-snp 50 --out results/qc/plink_roh
        python scripts/plot_roh.py results/qc/plink_roh.hom {output.plot}
        cp results/qc/plink_roh.hom {output.hom}
        rm results/qc/ann_roh.vcf.gz results/qc/temp_roh*
        """

rule genomic_windows:
    input: "results/qc/clean_full.vcf.gz"
    output: pi="results/qc/diversity.windowed.pi", tajima="results/qc/diversity.Tajima.D", plot="results/plots/Genomic_Windows.png"
    shell: "vcftools --gzvcf {input} --window-pi 50000 --out results/qc/diversity; vcftools --gzvcf {input} --TajimaD 50000 --out results/qc/diversity; python scripts/plot_windows.py {output.pi} {output.tajima} {output.plot}"
