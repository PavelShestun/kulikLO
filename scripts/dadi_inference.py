import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import dadi
import matplotlib.pyplot as plt
import numpy as np

vcf_file = sys.argv[1]
pop_map = sys.argv[2]
out_params = sys.argv[3]
out_plot = sys.argv[4]

# 1. Загрузка
dd = dadi.Misc.make_data_dict_vcf(vcf_file, pop_map)
# Используем проекцию поменьше (20), чтобы уменьшить шум отсутствующих данных
fs = dadi.Spectrum.from_data_dict(dd, ['Spoonbill'], [20], polarized=False)

# МАСКИРОВАНИЕ: Убираем синглетоны (ошибки секвенирования)
fs.mask[1] = True 
fs.mask[-1] = True

def bottleneck(params, ns, pts):
    nuB, nuF, TB, TF = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi, xx, TB, nuB)
    phi = dadi.Integration.one_pop(phi, xx, TF, nuF)
    return dadi.Spectrum.from_phi(phi, ns, (xx,))

pts_l = [30, 40, 50]
func_ex = dadi.Numerics.make_extrap_log_func(bottleneck)

p0 = [0.1, 1.0, 0.1, 0.1]
lower_bound = [1e-3, 1e-2, 1e-3, 1e-3]
upper_bound = [10, 100, 5, 5] # Расширили границы

best_ll = -np.inf
best_params = p0
best_model = None

# Быстрый прогон оптимизации
print("Optimizing...")
for i in range(15):
    p_guess = dadi.Misc.perturb_params(p0, fold=1)
    try:
        popt = dadi.Inference.optimize_log(p_guess, fs, func_ex, pts_l, 
                                           lower_bound=lower_bound, 
                                           upper_bound=upper_bound, 
                                           maxiter=20, verbose=0)
        model = func_ex(popt, fs.sample_sizes, pts_l)
        ll = dadi.Inference.ll_multinom(model, fs)
        if ll > best_ll:
            best_ll = ll
            best_params = popt
            best_model = model
    except: pass

if best_model is None:
    best_model = func_ex(p0, fs.sample_sizes, pts_l)

# РИСОВАНИЕ
fig = plt.figure(212, figsize=(8, 6))
dadi.Plotting.plot_1d_comp_multinom(best_model, fs)
plt.tight_layout()
plt.savefig(out_plot)

# Расчет (Theta based on unmasked sites)
theta = dadi.Inference.optimal_sfs_scaling(best_model, fs)
# N_anc calculation
mu = 4.6e-9
L = 5000000 # Примерная длина эффективного генома после фильтров
N_anc = theta / (4 * mu * L)
nuB, nuF, TB, TF = best_params

with open(out_params, 'w') as f:
    f.write(f"N_ANC={int(N_anc)}\n")
    f.write(f"N_BOT={int(N_anc*nuB)}\n")
    f.write(f"N_CUR={int(N_anc*nuF)}\n")
    f.write(f"T_BOT_DUR_GEN={int(TB*2*N_anc)}\n")
    f.write(f"T_REC_AGO_GEN={int(TF*2*N_anc)}\n")
