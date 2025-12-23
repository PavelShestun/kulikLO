import sys
import subprocess
import matplotlib.pyplot as plt

param_file = sys.argv[1]
out_plot = sys.argv[2]

# Чтение параметров
p = {}
with open(param_file) as f:
    for line in f:
        k, v = line.strip().split('=')
        p[k] = int(v)

# Scaling Factor (Q)
# Уменьшим Q, чтобы сохранить точность при малых N, но это замедлит расчет
Q = 5 

# --- ЛОГИКА ЗАЩИТЫ ОТ ВЫМИРАНИЯ ---
# Минимальный размер популяции в симуляции, чтобы она не вымерла случайно
MIN_VIABLE_SCALED = 20 

def scale_pop(val, q, min_val):
    scaled = int(val / q)
    return max(min_val, scaled)

N_ANC = scale_pop(p['N_ANC'], Q, MIN_VIABLE_SCALED)
# Здесь была проблема: 39/10 = 3 особи -> вымирание. Теперь будет минимум 20.
N_BOT = scale_pop(p['N_BOT'], Q, MIN_VIABLE_SCALED) 
N_CUR = scale_pop(p['N_CUR'], Q, MIN_VIABLE_SCALED)
N_CRISIS = scale_pop(2500, Q, MIN_VIABLE_SCALED)

# Время тоже скалируем, но следим, чтобы оно не стало 0
T_BOT = max(10, int(p['T_BOT_DUR_GEN'] / Q))
T_REC = max(10, int(p['T_REC_AGO_GEN'] / Q))

T_BURNIN = 10 * N_ANC 
T_START_BOT = T_BURNIN
T_START_REC = T_START_BOT + T_BOT
T_START_CRISIS = T_START_REC + T_REC
T_END = T_START_CRISIS + 100 

print(f"Simulating with scaled params: N_ANC={N_ANC}, N_BOT={N_BOT}, N_CUR={N_CUR}")

slim_code = f"""
initialize() {{
    initializeMutationRate(1e-7 * {Q}); 
    initializeMutationType("m1", 0.5, "f", 0.0);
    initializeMutationType("m2", 0.1, "g", -0.05, 0.2); 
    initializeGenomicElementType("g1", c(m1, m2), c(1.0, 0.5));
    for (i in 0:9) {{ initializeGenomicElement(g1, i*2000, i*2000 + 999); }}
    initializeRecombinationRate(1e-8 * {Q});
}}
1 early() {{ sim.addSubpop("p1", {N_ANC}); }}

1:{T_END} late() {{
    if (community.tick % 10 == 0) {{
        if (p1.individualCount > 0) {{
            load = mean(p1.individuals.countOfMutationsOfType(m2));
            cat("DATA:" + community.tick + "," + load + "\\n");
        }}
    }}
}}

{T_START_BOT} early() {{ p1.setSubpopulationSize({N_BOT}); }}
{T_START_REC} early() {{ p1.setSubpopulationSize({N_CUR}); }}
{T_START_CRISIS} early() {{ p1.setSubpopulationSize({N_CRISIS}); }}
"""

res = subprocess.run(["slim", "-"], input=slim_code, capture_output=True, text=True)

gens, loads = [], []
for line in res.stdout.splitlines():
    if line.startswith("DATA:"):
        try:
            parts = line.strip().split("DATA:")[1].split(",")
            gens.append(int(parts[0]))
            loads.append(float(parts[1]))
        except:
            continue

plt.figure(figsize=(10, 6))
if len(gens) > 0:
    plt.plot(gens, loads, 'r-', lw=2)
    plt.axvline(T_START_BOT, ls='--', color='black', label="Bottleneck")
    plt.axvline(T_START_REC, ls='--', color='blue', label="Recovery")
    plt.axvline(T_START_CRISIS, ls='--', color='red', label="Modern Crisis")
    plt.xlabel("Generations")
    plt.ylabel("Avg Load")
    plt.title(f"Simulated Load (Q={Q})")
    plt.legend()
else:
    plt.text(0.5, 0.5, "Simulation Extinct", ha='center')
    print("SLiM produced no data. Population likely extinct.")

plt.tight_layout()
plt.savefig(out_plot)
