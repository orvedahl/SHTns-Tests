"""
Read Fortran grid files and compare to other Python codes
"""
from __future__ import print_function
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.getenv("HOME"), "Programs", "RayleighUtils", "Utilities"))
from spectral_utils import Ylm

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

def compare(x, y):
    if (not np.allclose(x, y)):
        diff = x-y
        max_err = np.max(diff)
        min_err = np.min(diff)
        return False
    return True

def main():

    #---
    ylm = read("test_ylm.txt")
    N = np.shape(ylm)[0]
    results = []
    for i in range(N):
        l = int(ylm[i,0]); m = int(ylm[i,1])
        th = ylm[i,2]; phi = ylm[i,3]

        y = Ylm(l,m,th,phi)
        r = compare(np.real(y), ylm[i,4]); results.append(r)
        r = compare(np.imag(y), ylm[i,5]); results.append(r)

    if (not np.all(results)):
        print("\nFAILED\n")
    else:
        print("\nPASS\n")

if __name__ == "__main__":
    main()
