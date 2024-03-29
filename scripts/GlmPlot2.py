import numpy as np
import matplotlib.pyplot as plt
import os.path
import sys

from map_tools import *
from map_cycle import *
from map_run import *

#####################################################
### usage: python3 GlmPlot.py run_number1 run_number2
### saves a spectrum plot of the two runs in ../plots

run1 = int(sys.argv[1])
run2 = int(sys.argv[2])
ylim = int(sys.argv[3])
if len(sys.argv)>4:
    lmax0 = int(sys.argv[4])
else:
    lmax0 = 7

Gs1 = np.load("data/spectrum_run{}_lmax7.npy".format(run1), allow_pickle=True)
G1 = Gs1.item()['mainly_z']

Gs2 = np.load("data/spectrum_run{}_lmax7.npy".format(run2), allow_pickle=True)
G2 = Gs2.item()['mainly_z']

plt.rcParams["font.size"]=10
plt.rcParams["lines.markersize"]=3
plt.rcParams["figure.figsize"]=[12, 4*lmax0]
plt.rcParams["figure.dpi"]=200

dm=0.1
fig, axs = plt.subplots(lmax0, 1)
M = np.arange(-lmax0-1, lmax0+2)

for l in range(1, lmax0+1):
    l0 = l-1
    L = abs(vars(Geometry)['D{}'.format(l)])
    for m in M:
        if l0==0 and m==0:
            axs[l0].bar(m, 1e3*L*G1[l][7+1+m], color='green', alpha=0.5, width=dm)
            if len(sys.argv)>2:
                axs[l0].bar(m+dm, 1e3*L*G2[l][7+1+m], color='blue', alpha=0.5, width=dm)
        else:
            axs[l0].bar(m, 1e3*L*G1[l][7+1+m], color='green', alpha=0.5, width=dm)
            if len(sys.argv)>2:
                axs[l0].bar(m+dm, 1e3*L*G2[l][7+1+m], color='blue', alpha=0.5, width=dm)
    axs[l0].grid()
    axs[l0].set_axisbelow(True)
    #axs[l0].locator_params(axis='y', nbins=5)
    axs[l0].set_ylabel(r"$D_{{{}}}^{{{}}} \times G_{{ {}m }}$ (fT/cm)".format(l, l-1, l), size=20)
    axs[l0].set_ylim(-ylim, ylim)
    axs[l0].set_xticks(M)
    axs[l0].set_xticklabels([])
    if l==lmax0:
        axs[l0].set_xticklabels(M)
        axs[l0].set_xlabel(r"$m$", size=20)
axs[0].legend(fontsize=12)
axs[0].set_title(label=r"Runs {} (green) and {} (blue)".format(run1, run2), size=20)
# ax.text(0.1, 0.72, r'$I = -10$ mA', color='tab:blue', transform=ax.transAxes, size=12)
# ax.text(0.1, 0.8, r'$I = +10$ mA', color='tab:red', transform=ax.transAxes, size=12)
plt.savefig("plots/spectrumPlot_runs{}and{}_1.png".format(run1, run2), dpi=200)

# plt.show()
