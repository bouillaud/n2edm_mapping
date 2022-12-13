import numpy as np
import matplotlib.pyplot as plt
import os.path
import sys

from map_tools import *
from map_cycle import *
from map_run import *

################################################################################################
### usage: python3 GlmExtraction.py run_number maximum_l
### saves a numpy array containing Glm spectra from z, rho, and z-rho probes combined in ../data
### as well as a spectrum plot in ../plots

run = int(sys.argv[1])
lmax = int(sys.argv[2])

runmap = Run(run, np.arange(0, 150), dtformat='v4_PhiZswap', calib=True, trim=True)

### get Glms
G, G_err = runmap.getGlm(lmax, source='z')
G, G_err = runmap.getGlm(lmax, source='rho')
G, G_err = runmap.combineGlm(mode='z')

## save Glms
Gname = "spectrum_run{}_lmax{}".format(run, lmax)
Gerrname = "spectrumErr_run{}_lmax{}".format(run, lmax)
#i=0
#while os.path.exists("../data/{}_{}.npy".format(Gname, i)):
    #i+=1
#np.save("../data/{}_{}".format(Gname, i), runmap.Gs, allow_pickle=True)
#np.save("../data/{}_{}".format(Gerrname, i), runmap.Gs_err, allow_pickle=True)
np.save("../data/{}".format(Gname), runmap.Gs, allow_pickle=True)
np.save("../data/{}".format(Gerrname), runmap.Gs_err, allow_pickle=True)
print("Gs saved")

print("G_00 = {:.2e} +/- {:.2e}".format(G[0][lmax+1+0], G_err[0][lmax+1+0]))
print("G_1 = {:.2e} +/- {:.2e}".format(G[1][lmax+1+0], G_err[1][lmax+1+0]))
print("phanG_3 = {:.2e} +/- {:.2e}".format(G[3][lmax+1+0]*Geometry.L3/Geometry.C3, G_err[3][lmax+1+0]**Geometry.L3/Geometry.C3))
print("phanG_5 = {:.2e} +/- {:.2e}".format(G[5][lmax+1+0]*Geometry.L5/Geometry.C5, G_err[5][lmax+1+0]*Geometry.L5/Geometry.C5))

### plot Glms

dm=0.1
lmax0=lmax
plt.rcParams["font.size"]=10
plt.rcParams["lines.markersize"]=3
plt.rcParams["figure.figsize"]=[12, 4*lmax0]
plt.rcParams["figure.dpi"]=200

fig, axs = plt.subplots(lmax0, 1)
M = np.arange(-lmax0-1, lmax0+2)

for l in range(1, lmax0+1):
    l0 = l-1
    L = abs(vars(Geometry)['D{}'.format(l)])
    #print(L)
    for m in M:
        if l0==0 and m==0:
            axs[l0].bar(m, L*G[l][7+1+m], yerr=L*G_err[l][7+1+m], align='center', color='tab:blue', ecolor='black', capsize=5, width=dm)
        else:
            axs[l0].bar(m, L*G[l][7+1+m], yerr=L*G_err[l][7+1+m], align='center', color='tab:blue', ecolor='black', capsize=5, width=dm)
    axs[l0].grid()
    #axs[l0].locator_params(axis='y', nbins=5)
    axs[l0].set_ylabel(r"$D_{{{}}}^{{{}}} \times G_{{ {}m }}$ (pT/cm)".format(l, l-1, l), size=20)
    #axs[l0].set_ylim(-4, 4)
    axs[l0].set_xticks(M)
    axs[l0].set_xticklabels([])
    if l==lmax0:
        axs[l0].set_xticklabels(M)
        axs[l0].set_xlabel(r"$m$", size=20)
axs[0].legend(fontsize=12)
axs[0].set_title(label=r"Run {}".format(run), size=20)
# ax.text(0.1, 0.72, r'$I = -10$ mA', color='tab:blue', transform=ax.transAxes, size=12)
# ax.text(0.1, 0.8, r'$I = +10$ mA', color='tab:red', transform=ax.transAxes, size=12)
plt.savefig("../plots/spectrumPlot_run{}_lmax{}.png".format(run, lmax0), dpi=200)
