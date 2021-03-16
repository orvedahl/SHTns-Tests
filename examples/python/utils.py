"""
Various useful routines for interfacing with Rayleigh
"""

from __future__ import print_function
import numpy as np
from scipy.special import factorial as Factorial
import NumericalTools.public as NT

def Ylm(ell, m, theta, phi):
    """
    explicitly evaluate the spherical harmonic on the given (theta,phi) grid,
    does not rely on the recurrence formula

    either theta or phi can be arrays, just not both
    """
    if (m == 0): # catch the m=0 cases
        return Yl0(ell, theta)

    if (abs(m) > ell):
        raise ValueError("Ylm must have |m| <= l; l={}, m={}".format(ell,m))

    abs_m = abs(m)
    ell_m = 1000000*ell + abs_m

    x = np.cos(theta)
    y = np.sin(theta) # y = (1-x**2)**0.5

    # m/=0 and m>0 cases
    if (ell_m == 1000001): # l=1,m=1
        Plm = -y
    elif (ell_m == 2000001): # l=2,m=1
        Plm = -3*x*y
    elif (ell_m == 2000002): # l=2,m=2
        Plm = 3*y*y
    elif (ell_m == 3000001): # l=3,m=1
        Plm = -1.5*(5*x*x-1.)*y
    elif (ell_m == 3000002): # l=3,m=2
        Plm = 15*x*y*y
    elif (ell_m == 3000003): # l=3,m=3
        Plm = -15*y*y*y
    elif (ell_m == 4000001): # l=4,m=1
        Plm = -2.5*x*(7.*x*x - 3.)*y
    elif (ell_m == 4000002): # l=4,m=2
        Plm = 7.5*(7.*x*x - 1.)*y*y
    elif (ell_m == 4000003): # l=4,m=3
        Plm = -105*x*y*y*y
    elif (ell_m == 4000004): # l=4,m=4
        Plm = 105*y*y*y*y
    else:
        raise ValueError("Ylm has not been coded for \ell={} and m>0".format(ell))

    ylm = Alm(ell,abs_m)*Plm*np.exp(1j*abs_m*phi)

    # catch negative m values
    if (m < 0):
        ylm = (-1)**m * np.conjugate(ylm)

    return ylm

def Yl0(ell, theta):
    """
    explicityl compute some m=0 spherical harmonics, does not rely on the recurrence formula

    see https://en.wikipedia.org/wiki/Table_of_spherical_harmonics
    """
    x = np.cos(theta)
    if (ell == 0):
        Pl = 1.0
    elif (ell == 1):
        Pl = x
    elif (ell == 2):
        Pl = 0.5*(3*x*x - 1)
    elif (ell == 3):
        Pl = 0.5*(5*x**3 - 3*x)
    elif (ell == 4):
        Pl = (35*x**4 - 30*x**2 + 3.)/8.
    elif (ell == 5):
        Pl = (63*x**5 - 70*x**3 + 15*x)/8.
    elif (ell == 6):
        Pl = (231*x**6 - 315*x**4 + 105*x**2 - 5)/16.
    elif (ell == 7):
        Pl = (429*x**7 - 693*x**5 + 315*x**3 - 35*x)/16.
    elif (ell == 8):
        Pl = (6435*x**8 - 12012*x**6 + 6930*x**4 - 1260*x**2 + 35)/128.
    elif (ell == 9):
        Pl = (12155*x**9 - 25740*x**7 + 18018*x**5 - 4620*x**3 + 315*x)/128.
    elif (ell == 10):
        Pl = (46189*x**10 - 109395*x**8 + 90090*x**6 - 30030*x**4 + 3465*x**2 - 63)/256.
    else:
        raise ValueError("Yl0 has not been coded for \ell={}".format(ell))
    return Alm(ell,0)*Pl

def Alm(l,m):
    """
    spherical harmonic normalization
       Alm**2 = (2*l+1)/4/pi * (l-m)! / (l+m)!
    l is assumed >0
    m can be negative
    """
    A = (2*l+1.)/4/np.pi
    A *= Factorial(l-m)
    A /= Factorial(l+m)
    A = np.sqrt(A)
    return A

class LegendreTransform():
    """
    streamline the transform process: AzAvg <--> m=0 spectra
    """
    def __init__(self, N, spectral, dealias=1.5):
        """
        N - integer, either nth or lmax
        spectral - bool, is N=nth or N=lmax
        """
        if (spectral):
            lmax = N
            nth = int(dealias*(lmax+1))
        else:
            nth = N
            lmax = int(nth/dealias - 1)

        self.LT = NT.Transforms.LegendreTransform(nth, a=-1, b=1, dealias=dealias)
        self.nth = self.LT.Ngrid
        self.nl = self.LT.Npoly
        self.lmax = self.LT.Npoly_max

    def ToSpectral(self, data, axis=0):
        return self.LT.ToSpectral(data, axis=axis)

    def ToPhysical(self, data, axis=0):
        return self.LT.ToPhysical(data, axis=axis)

    def Finalize(self):
        self.Close()
    def Close(self):
        self.LT.Finalize()

def RadialDomain(L, chi):
    """
    calculate radial domain based on shell depth and aspect ratio
    """
    r_outer = L/(1.-chi)
    r_inner = r_outer*chi
    return r_inner, r_outer

def RadialGrid(nr, domain, physical_coords=False):
    """
    build radial Rayleigh grid

    nr     --- number of radial grid points
    domain --- bounds of radial grid

    physical_coords --- if True, "domain" is a list holding [rin, rout]
                        if False, "domain" holds [chi, L]
    """
    # angular parameters do not mater, so just pick a small number
    R = RayleighGrid(nr, 16, domain, send_lmax=False, physical_coords=physical_coords)
    radius = R.radius
    return radius

def AngularGrid(n, send_lmax=True, degrees=False):
    """
    build angular Rayleigh grid

    n --- resolution of angular grid

    send_lmax --- if True, "n" refers to the maximum spherical harmonic
                  if False, "n" is the number of physical grid points
    """
    # radial parameters do not mater, so just pick a small number
    R = RayleighGrid(16, n, [0.35,1.], send_lmax=send_lmax, physical_coords=False)
    if (degrees):
        th  = R.theta_deg
        phi = R.phi_deg
    else:
        th  = R.theta
        phi = R.phi
    return th, phi

class RayleighGrid():
    """
    Attributes
    ----------
    lmax, nth, nphi, nr - integers specifying grid size
    r_inner, r_outer - floats specifying radial domain bounds
    radius - radial grid
    costheta, sintheta - grids in latitude
    theta - latitude grid (radians)
    phi - longitude grid (radians)
    """

    def __init__(self, nr, n, domain, send_lmax=True, physical_coords=False):
        """
        build a Rayleigh grid

        nr     --- number of radial grid points
        n      --- resolution of angular grid
        domain --- bounds of radial grid

        send_lmax       --- if True, "n" refers to the maximum spherical harmonic
                            if False, "n" is the number of physical grid points

        physical_coords --- if True, "domain" is a list holding [rin, rout]
                            if False, "domain" holds [chi, L]
        """

        if (send_lmax):
            # incoming "n" is lmax, calculate nth
            self.lmax = int(n)
            self.nth  = int(3./2.*(self.lmax+1))
        else:
            # incoming "n" is nth, calculate lmax
            self.nth  = int(n)
            self.lmax = int(2./3.*self.nth - 1)
        self.nphi = 2*self.nth
        self.nr = int(nr)

        # get bounds of radial grid
        if (physical_coords):
            # physical coordinates were given
            r_inner = domain[0]
            r_outer = domain[1]
        else:
            # nondimensional values were given
            chi = domain[0]
            L   = domain[1]
            r_inner, r_outer = RadialDomain(L, chi)
        self.r_inner = r_inner
        self.r_outer = r_outer
        self.L = self.r_outer - self.r_inner
        self.chi = self.r_inner / self.r_outer

        self.domain = [[self.r_inner, self.r_outer], [0, np.pi], [0, 2*np.pi]]

        # build/extract radial grid
        self.radial_grid = NT.ChebGrid(self.nr, a=self.r_inner, b=self.r_outer,
                                       zeros=True, endpoints=True)
        # reverse radial grid
        self.radius = self.radial_grid.xab[::-1]

        # build/extract angular grid
        self.costheta_grid = NT.LegendreGrid(self.nth, a=-1.0, b=1.0, endpoints=False)
        self.costheta = self.costheta_grid.xab   # array increases from -1 to 1
        self.sintheta = 1. - self.costheta**2
        self.theta    = np.arccos(self.costheta) # array increases from pi to 0
        self.latitude = np.pi/2. - self.theta
        self.theta_deg    = self.theta*180./np.pi
        self.latitude_deg = self.latitude*180./np.pi

        self.phi_grid = NT.UniformGrid(self.nphi, a=0, b=2*np.pi,
                                       lendpoint=True, uendpoint=False)
        self.phi       = self.phi_grid.x
        self.longitude = self.phi
        self.phi_deg       = self.phi*180./np.pi
        self.longitude_deg = self.longitude*180./np.pi

