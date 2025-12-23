args <- commandArgs(trailingOnly = TRUE)
library(LEA)

input_vcf <- args[1]
k_min     <- as.integer(args[2])
k_max     <- as.integer(args[3])
reps      <- as.integer(args[4])
out_plot  <- args[5]

# Конвертация
geno_file <- vcf2geno(input_vcf, output = paste0(input_vcf, ".geno"), force=TRUE)

# Запуск
project <- snmf(geno_file, K = k_min:k_max, entropy = TRUE, repetitions = reps, project = "new")

# Просто рисуем график энтропии (без ручного поиска минимума, на котором был сбой)
png(out_plot, width=800, height=600)
plot(project, col = "blue", pch = 19, type = "b")
dev.off()
