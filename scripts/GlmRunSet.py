import numpy as np
import matplotlib.pyplot as plt
import os.path
import sys

from map_tools import *
from map_cycle import *
from map_run import *

################################################################################################
### usage: python3 GlmRunSet.py run_start run_stop maximum_l dtformat

runstart = int(sys.argv[1])
runstop = int(sys.argv[2])
lmax = int(sys.argv[3])
dtformat = str(sys.argv[4])

G10s = []
for run in range(runstart, runstop):
    runmap = Run(run, np.arange(0, 150), dtformat=dtformat, calib=True, trim=True)
    ### get Glms
    G, G_err = runmap.getGlm(lmax, source='z')
    G10s.append(G[1][lmax+1])
print(G10s)

filename = "G10s_runstart{}_runstop{}_lmax{}".format(runstart, runstop, lmax)
np.save("../data/{}".format(filename), G10s, allow_pickle=True)


#G, G_err = runmap.getGlm(lmax, source='rho')
#G, G_err = runmap.combineGlm(mode='z')

## save Glms

#Gname = "spectrum_run{}_lmax{}".format(run, lmax)
#Gerrname = "spectrumErr_run{}_lmax{}".format(run, lmax)
##i=0
##while os.path.exists("../data/{}_{}.npy".format(Gname, i)):
#    #i+=1
##np.save("../data/{}_{}".format(Gname, i), runmap.Gs, allow_pickle=True)
##np.save("../data/{}_{}".format(Gerrname, i), runmap.Gs_err, allow_pickle=True)
#np.save("../data/{}".format(Gname), runmap.Gs, allow_pickle=True)
#np.save("../data/{}".format(Gerrname), runmap.Gs_err, allow_pickle=True)
#print("Gs saved")
#print("G_00 = {:.2e} +/- {:.2e}".format(G[0][lmax+1+0], G_err[0][lmax+1+0]))
#print("G_1 = {:.2e} +/- {:.2e}".format(G[1][lmax+1+0], G_err[1][lmax+1+0]))
#if lmax>=3:
#    print("phanG_3 = {:.2e} +/- {:.2e}".format(G[3][lmax+1+0]*Geometry.L3/Geometry.C3, G_err[3][lmax+1+0]**Geometry.L3/Geometry.C3))
#if lmax>=5:
#    print("phanG_5 = {:.2e} +/- {:.2e}".format(G[5][lmax+1+0]*Geometry.L5/Geometry.C5, G_err[5][lmax+1+0]*Geometry.L5/Geometry.C5))
#if lmax>=7:
#    print("phanG_7 = {:.2e} +/- {:.2e}".format(G[7][lmax+1+0]*Geometry.L7/Geometry.C7, G_err[7][lmax+1+0]*Geometry.L7/Geometry.C7))
#
