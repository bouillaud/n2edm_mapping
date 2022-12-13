import numpy as np
import matplotlib.pyplot as plt
import os.path
import sys

from map_tools import *
from map_cycle import *
from map_run import *


### edit plot parameters here:
plt.rcParams["font.size"]=10
plt.rcParams["lines.markersize"]=3
plt.rcParams["figure.figsize"]=[12, 6]
plt.rcParams["figure.dpi"]=100

############################
### 3 SCRIPT ARGUMENTS
run = int(sys.argv[1])
rho = int(sys.argv[2])
z = int(sys.argv[3])
Bi = str(sys.argv[4]) # field projection Brho Bphi or Bz

runmap = Run(run, np.arange(0, 150), dtformat='v4_PhiZswap', calib=True, trim=True)

ring = runmap.cycles['rings'][rho, z]
print("This is a {} scan at rho={:.0f}cm, z={:.0f}cm".format(ring.scantype, ring.rho[0], ring.z[0]))

fig, ax = ring.simplePlot('phi', Bi, fontsize=15)
plt.savefig("../plots/ringplot_run{}_rho{}_z{}_{}.png".format(run, rho, z, Bi), dpi=200)

plt.show()
