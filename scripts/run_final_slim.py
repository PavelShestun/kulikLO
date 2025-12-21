import os
import subprocess
import matplotlib.pyplot as plt
import sys

OUT_PLOT = sys.argv[1]

# Параметры (как в вашем примере)
N_ANC = 68440; N_BOT = 20; N_CUR = 16530
GEN_TIME = 4
Q = 10 # Scaling

N_ANC_S = int(N_ANC / Q); N_BOT_S = int(N_BOT); N_CUR_S = int(N_CUR / Q)
T_BOT_DUR_S = int((10895 / GEN_TIME) / Q)
T_REC_AGO_S = int((114219 / GEN_TIME) / Q)
T_START_BOT = 5000
T_START_REC = T_START_BOT + T_BOT_DUR_S
T_END = T_START_REC + T_REC_AGO_S

slim_code = f"""
initialize() {{
    initializeMutationRate(1e-7 * {Q}); 
    initializeMutationType("m1", 0.5, "f", 0.0);
    initializeMutationType("m2", 0.05, "g", -0.05, 0.2); 
    initializeGenomicElementType("g1", c(m1, m2), c(1.0, 2.0));
    for (i in 0:49) {{ initializeGenomicElement(g1, i*2000, i*2000 + 999); }}
    initializeRecombinationRate(1e-8 * {Q});
}}
1 early() {{ sim.addSubpop("p1", {N_ANC_S}); }}
1:{T_END} late() {{
    if (community.tick % 100 == 0) {{
        if (p1.individualCount > 0) {{
            inds = p1.individuals;
            total_m2 = mean(inds.countOfMutationsOfType(m2));
            cat(community.tick + "," + total_m2 + "\\n");
        }}
    }}
}}
{T_START_BOT} early() {{ p1.setSubpopulationSize({N_BOT_S}); }}
{T_START_REC} early() {{ p1.setSubpopulationSize({N_CUR_S}); }}
"""

with open("temp_model.slim", "w") as f: f.write(slim_code)
res = subprocess.run(["slim", "temp_model.slim"], capture_output=True, text=True)

gens, loads = [], []
for line in res.stdout.splitlines():
    if "," in line and not line.startswith("Gen"):
        p = line.split(",")
        try: gens.append(int(p[0])); loads.append(float(p[1]))
        except: pass

if len(gens) > 0:
    plt.figure(figsize=(10, 6))
    plt.plot(gens, loads, color="red")
    plt.axvline(x=T_START_BOT, linestyle='--'); plt.axvline(x=T_START_REC, linestyle='--')
    plt.savefig(OUT_PLOT)
