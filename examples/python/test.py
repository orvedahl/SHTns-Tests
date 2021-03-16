from __future__ import print_function
import utils
import numpy as np
import NumericalTools.public as NT

def ToSpectral(fx, iPlm):
    """
    fx = (nx,)
    iPlm = (nx,nl), appropriate for chosen m

    Clm = (nl,)
    """
    return np.dot(np.transpose(iPlm), fx)

def ToPhysical(Clm, Plm):
    """
    Clm = (nl,)
    Plm = (nx,nl), appropriate for chosen m

    fx = (nx,)
    """
    return np.dot(Plm,Clm)

def main():
    nr = 16
    lmax = 31

    Rgrid = utils.RayleighGrid(nr, lmax, [0.35, 1.0])
    costh = Rgrid.costheta
    sinth = Rgrid.sintheta
    theta = Rgrid.theta
    phi = Rgrid.phi

    weights = Rgrid.costheta_grid.weights
    Weights = np.reshape(weights, np.shape(weights)+(1,1))

    Plm = NT.AlmPlm(lmax, True) # (x,l,m)
    iPlm = 2*np.pi*Weights*Plm

    l = 3; m = 2
    true_phys = np.real(utils.Ylm(l, m, theta, 0.0))
    true_spec = np.zeros((lmax+1,))
    true_spec[l] = 1.0

    spectral = ToSpectral(true_phys, iPlm[:,:,m])
    physical = ToPhysical(spectral, Plm[:,:,m])

    err = np.max(np.abs(physical - true_phys))
    print("m={}, max err in Phys-->Spec-->Phys: {}".format(m,err))

    maxerr = np.max(np.abs(spectral - true_spec))
    print("m={}, max err in Phys-->Spec: {}".format(m,maxerr))

    physical = ToPhysical(true_spec, Plm[:,:,m])
    spectral = ToSpectral(physical, iPlm[:,:,m])

    err = np.max(np.abs(spectral - true_spec))
    print("m={}, max err in Spec-->Phys-->Spec: {}".format(m,err))

    maxerr = np.max(np.abs(physical - true_phys))
    print("m={}, max err in Spec-->Phys: {}".format(m,maxerr))

if __name__ == "__main__":
    main()
