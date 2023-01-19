import numpy as np
import matplotlib.pyplot as plt
import os.path
import sys

### takes 2 arguments: run number and field projection ###
nr = int(sys.argv[1])
iB = int(sys.argv[2])

########################################
### edit plot parameters here
plt.rcParams["font.size"]=10
plt.rcParams["lines.markersize"]=3
plt.rcParams["figure.figsize"]=[12, 12]
plt.rcParams["figure.dpi"]=200
########################################

def runPlot(Rho, Phi, Z, Brho, BPhi, Bz, iB=2, cmap='magma', zratio=6, elev=15, azim=30, clim=[]):
    """
    takes: 3 2D arrays of scan coordinates and 3 2D arrays of field projections, in cyclindrical coordinates, and index iB=0,1,2 which selects the field projection to plot.
    returns: figure and ax object of a 3D plot.
    """

    X = Rho*np.cos(Phi*np.pi/180)
    Y = Rho*np.sin(Phi*np.pi/180)

    fig = plt.figure(figsize=plt.figaspect(1))
    ax = fig.add_subplot(projection='3d')

    if iB==0:
        p = ax.scatter(X, Y, Z, c=Brho, cmap=cmap)
        cb = fig.colorbar(p, ax=ax, shrink=0.5)
        cb.set_label(r"$B_\rho$ (pT)", size=15)
    elif iB==1:
        p = ax.scatter(X, Y, Z, c=Bphi, cmap=cmap)
        cb = fig.colorbar(p, ax=ax, shrink=0.5)
        cb.set_label(r"$B_\varphi$ (pT)", size=15)
    elif iB==2:
        p = ax.scatter(X, Y, Z, c=Bz, cmap=cmap)
        cb = fig.colorbar(p, ax=ax, shrink=0.5)
        cb.set_label(r"$B_z$ (pT)", size=15)

    ax.set_xlabel(r"$x$ (cm)", size=12)
    ax.set_ylabel(r"$y$ (cm)", size=12)
    ax.set_zlabel(r"$z$ (cm)", size=12)

    ax.view_init(elev, azim)
    ax.set_box_aspect((4, 4, zratio))
    if len(clim)>1:
        p.set_clim(clim[0], clim[1])

    return fig, ax


Rho = np.loadtxt('data/B0map{}_rho.txt'.format(nr))
Phi = np.loadtxt('data/B0map{}_phi.txt'.format(nr))
Z = np.loadtxt('data/B0map{}_z.txt'.format(nr))
Brho = np.loadtxt('data/B0map{}_Brho.txt'.format(nr))
Bphi = np.loadtxt('data/B0map{}_Bphi.txt'.format(nr))
Bz = np.loadtxt('data/B0map{}_Bz.txt'.format(nr))

runPlot(Rho, Phi, Z, Brho, Bphi, Bz, iB)
plt.savefig('runPlot_run{}_iB{}.png'.format(nr, iB))
