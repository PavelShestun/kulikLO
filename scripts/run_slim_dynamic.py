import sys
import subprocess
import matplotlib.pyplot as plt

param_file = sys.argv[1]
out_plot = sys.argv[2]

p = {}
with open(param_file) as f:
    for line in f:
        k, v = line.strip().split('=')
        p[k] = int(v)

Q = 500 # Оптимальное ускорение
N_ANC = int(p['N_ANC']/Q); N_BOT = int(p['N_BOT']/Q); N_CUR = int(p['N_CUR']/Q)
T_BOT = int(p['T_BOT_DUR_GEN']/Q); T_REC = int(p['T_REC_AGO_GEN']/Q)

# Современный кризис (реальные 25 000 птиц -> Ne примерно 2500)
N_CRISIS = int(2500 / Q)
if N_CRISIS < 2: N_CRISIS = 2

T_BURNIN = 500
T_START_BOT = T_BURNIN
T_START_REC = T_START_BOT + T_BOT
T_START_CRISIS = T_START_REC + T_REC
T_END = T_START_CRISIS + 100 # Моделируем 100 поколений кризиса

slim_code = f"""
initialize() {{
    initializeMutationRate(1e-7 * {Q}); 
    initializeMutationType("m1", 0.5, "f", 0.0);
    initializeMutationType("m2", 0.1, "g", -0.05, 0.2); // Вредные
    initializeGenomicElementType("g1", c(m1, m2), c(1.0, 0.5));
    for (i in 0:9) {{ initializeGenomicElement(g1, i*2000, i*2000 + 999); }}
    initializeRecombinationRate(1e-8 * {Q});
}}
1 early() {{ sim.addSubpop("p1", {N_ANC}); }}
1:{T_END} late() {{
    if (community.tick % 10 == 0) {{
        load = mean(p1.individuals.countOfMutationsOfType(m2));
        cat(community.tick + "," + load + "\\n");
    }}
}}
{T_START_BOT} early() {{ p1.setSubpopulationSize({N_BOT}); }}
{T_START_REC} early() {{ p1.setSubpopulationSize({N_CUR}); }}
{T_START_CRISIS} early() {{ p1.setSubpopulationSize({N_CRISIS}); }}
"""

with open("model.slim", "w") as f: f.write(slim_code)
res = subprocess.run(["slim", "model.slim"], capture_output=True, text=True)

gens, loads = [], []
for line in res.stdout.splitlines():
    if "," in line:
        parts = line.split(","); gens.append(int(parts[0])); loads.append(float(parts[1]))

plt.figure(figsize=(10, 6))
plt.plot(gens, loads, 'r-', lw=2)
plt.axvline(T_START_BOT, ls='--', color='black', label="Ancient Bottleneck")
plt.axvline(T_START_REC, ls='--', color='blue', label="Recovery")
plt.axvline(T_START_CRISIS, ls='--', color='red', label="Modern Crisis")
plt.title("Evolution of Mutational Load (with Modern Crisis)")
plt.legend(); plt.savefig(out_plot)
