import sys
import dadi
import matplotlib.pyplot as plt
import numpy as np

vcf_file = sys.argv[1]
pop_map = sys.argv[2]
out_params = sys.argv[3]
out_plot = sys.argv[4]

# Загрузка данных
dd = dadi.Misc.make_data_dict_vcf(vcf_file, pop_map)
# Проекция на 40 гаплотипов (для 22 особей)
fs = dadi.Spectrum.from_data_dict(dd, ['Spoonbill'], [40], polarized=False)

def bottleneck(params, ns, pts):
    nuB, nuF, TB, TF = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi, xx, TB, nuB)
    phi = dadi.Integration.one_pop(phi, xx, TF, nuF)
    return dadi.Spectrum.from_phi(phi, ns, (xx,))

pts_l = [40, 50, 60]
func_ex = dadi.Numerics.make_extrap_log_func(bottleneck)
# Запуск оптимизации
popt = dadi.Inference.optimize_log([0.1, 1.1, 0.1, 0.1], fs, func_ex, pts_l, 
                                   lower_bound=[1e-3, 1e-2, 1e-3, 1e-3], upper_bound=[1, 10, 1, 1], maxiter=30)

# РИСОВАНИЕ (Исправлено)
model = func_ex(popt, fs.sample_sizes, pts_l)
fig = plt.figure(figsize=(10, 8))
dadi.Plotting.plot_1d_comp_multinom(model, fs)
plt.savefig(out_plot, dpi=150)

# Расчет параметров
theta = dadi.Inference.optimal_sfs_scaling(model, fs)
N_anc = theta / (4 * 4.6e-9 * 1e8)
nuB, nuF, TB, TF = popt

with open(out_params, 'w') as f:
    f.write(f"N_ANC={int(N_anc)}\nN_BOT={int(N_anc*nuB)}\nN_CUR={int(N_anc*nuF)}\n")
    f.write(f"T_BOT_DUR_GEN={int(TB*2*N_anc)}\nT_REC_AGO_GEN={int(TF*2*N_anc)}\n")
