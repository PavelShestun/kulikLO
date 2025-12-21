args <- commandArgs(trailingOnly = TRUE)
input_vcf <- args[1]
k_min <- as.integer(args[2])
k_max <- as.integer(args[3])
reps <- as.integer(args[4])
out_plot <- args[5]

if (!requireNamespace("LEA", quietly = TRUE)) {
    stop("Package LEA is not installed. Run: conda install -c bioconda bioconductor-lea")
}
library(LEA)

# Конвертация VCF в формат GENO (требуется для sNMF)
# output будет input_vcf.geno
geno_file <- vcf2geno(input_vcf, output = paste0(input_vcf, ".geno"), force=TRUE)

# Запуск sNMF
project <- snmf(geno_file, K = k_min:k_max, entropy = TRUE, repetitions = reps, project = "new")

# Рисуем график кросс-энтропии
png(out_plot)
plot(project, col = "blue", pch = 19, cex = 1.2)
dev.off()

# Лучший K (минимальная энтропия)
best_K <- which.min(cross.entropy(project, K = k_min:k_max))
print(paste("Best K seems to be:", best_K))
