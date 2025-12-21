import os

# Имя итогового файла
output_filename = "full_pipeline_code.txt"

# Список файлов и папок, которые мы хотим собрать
# Мы берем только текстовые файлы с кодом и настройками
targets = [
    "Snakefile",
    "config/config.yaml",
    "environment.yaml",
    "scripts",  # Вся папка со скриптами
    "README.md"
]

# Расширения файлов, которые мы считаем "кодом"
code_extensions = (".py", ".R", ".slim", ".yaml", ".md", "") # "" для Snakefile

with open(output_filename, "w", encoding="utf-8") as out_file:
    for item in targets:
        if not os.path.exists(item):
            continue

        # Если это файл (например, Snakefile)
        if os.path.isfile(item):
            out_file.write(f"{item}:\n")
            with open(item, "r", encoding="utf-8") as f:
                out_file.write(f.read())
            out_file.write("\n\n" + "="*50 + "\n\n")

        # Если это папка (scripts/)
        elif os.path.isdir(item):
            for root, dirs, files in os.walk(item):
                for file in files:
                    if file.endswith(code_extensions):
                        full_path = os.path.join(root, file)
                        out_file.write(f"{full_path}:\n")
                        try:
                            with open(full_path, "r", encoding="utf-8") as f:
                                out_file.write(f.read())
                        except Exception as e:
                            out_file.write(f"Error reading file: {e}")
                        out_file.write("\n\n" + "="*50 + "\n\n")

print(f"Готово! Весь код собран в файл: {output_filename}")
