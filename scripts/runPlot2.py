import numpy as np
import matplotlib.pyplot as plt
import os.path
import sys

from map_tools import *
from map_cycle import *
from map_run import *

### arguments #######
nr = int(sys.argv[1]) # run number
dtformat = str(sys.argv[2]) # run data format
iB = int(sys.argv[3]) # field projection
cmap = str(sys.argv[4]) # plt colormap
ms = int(sys.argv[5]) # markersize

########################################
### edit plot parameters here
plt.rcParams["font.size"]=10
plt.rcParams["lines.markersize"]=3
plt.rcParams["figure.figsize"]=[12, 12]
plt.rcParams["figure.dpi"]=200
########################################

run = Run(nr, np.arange(0, 150), dtformat=dtformat, trim=True, calib=True)

run.newRunPlot(iB=iB, ms=ms, cmap=cmap)
plt.savefig('./plots/runPlot_nr{}_iB_{}'.format(nr, iB), dpi=200)
plt.show()
