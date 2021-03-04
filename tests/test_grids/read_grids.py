"""
Read Fortran grid files and compare to other Python codes
"""
from __future__ import print_function
import NumericalTools.public as NT
import numpy as np

def read(fname):
    entries = []
    with open(fname, "r") as f:
        for line in f:
            if (not line.lstrip().startswith("#")):
                fields = line.split()
                fields = [float(f) for f in fields]
                entries.append(fields)
    entries = np.array(entries) # (Nx, Nrows)
    return entries

def compare(x, y, label=""):
    if (len(x) != len(y)):
        print("Values are different sizes, {}".format(label))
        return
    if (not np.allclose(x, y)):
        print("Values are different, {}".format(label))
        diff = x-y
        max_err = np.max(diff)
        min_err = np.min(diff)
        print("\tmin error = {}".format(min_err))
        print("\tmax error = {}".format(min_err))
        return
    print("Values are the same, {}".format(label))

def main():

    #---
    cheb = read("test_grid_chebyshev.txt")
    N = np.shape(cheb)[0]
    Cheb = NT.ChebGrid(N, a=-1, b=1, zeros=1, endpoints=1) # Rayleigh grid is reversed

    compare(cheb[:,1], Cheb.x[::-1], "cheb x")
    compare(cheb[:,2], Cheb.theta[::-1], "cheb th")

    #---
    phi = read("test_grid_fourier.txt")
    N = np.shape(phi)[0]
    Phi = NT.UniformGrid(N, a=0, b=2*np.pi, lendpoint=1, uendpoint=0)

    compare(phi[:,1], Phi.x, "fourier")

    #---
    leg = read("test_grid_legendre.txt")
    N = np.shape(leg)[0]
    Leg = NT.LegendreGrid(N, a=-1, b=1, endpoints=0)
    x = Leg.xab
    th = np.arccos(x)
    w = Leg.weights

    compare(leg[:,1], x, "legendre x")
    compare(leg[:,2], th, "legendre th")
    compare(leg[:,3], w, "legendre weights")

    #---
    plm = read("test_grid_Plm.txt")
    N = np.shape(plm)[0]
    Plm = NT.AlmPlm(N, spectral=False) # (x, l, m)

    compare(plm[:,1], Plm[:,10,5], "Plm(:,10,5)")
    compare(plm[:,2], Plm[:,5,10], "Plm(:,5,10)")
    compare(plm[:,3], Plm[:,1,0], "Plm(:,1,0)")
    compare(plm[:,4], Plm[:,1,1], "Plm(:,1,1)")

if __name__ == "__main__":
    main()
