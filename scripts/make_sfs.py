import dadi
import matplotlib.pyplot as plt
import sys

vcf_file = sys.argv[1]
pop_file = sys.argv[2]
out_sfs = sys.argv[3]
out_plot = sys.argv[4]

dd = dadi.Misc.make_data_dict_vcf(vcf_file, pop_file)
pop_ids = ['Spoonbill']
proj_n = [40] 

fs = dadi.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections=proj_n, polarized=False)
fs.to_file(out_sfs)

plt.figure(figsize=(10, 6))
dadi.Plotting.plot_1d_fs(fs)
plt.title("SFS Spoonbill")
plt.savefig(out_plot)
