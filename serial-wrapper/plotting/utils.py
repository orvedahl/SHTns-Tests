from __future__ import print_function
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.rcParams.update({"text.usetex":True})
from matplotlib import cbook
from matplotlib.colors import Normalize, ListedColormap
import numpy as np
from numpy import ma

def ColorMap(cmap):
    """
    Build a colormap
    """
    if (cmap == 'default'):
        return mpl.cm.get_cmap('RdBu_r')
    elif (cmap == 'custom'):
        top = mpl.cm.get_cmap('Reds_r')
        bot = mpl.cm.get_cmap('Blues')

        N = 10
        tlo = 0.5; thi = 0.9
        blo = 0.5; bhi = 0.9

        newcols = np.vstack((top(np.linspace(tlo,thi,N)), bot(np.linspace(blo,bhi,N))))

        newcmp = ListedColormap(newcols, name='Custom')

        return newcmp

    else:
        raise ValueError("Unrecognized colormap choice: \"{}\"".format(cmap))

def MakeUnique(data, il, iNr, iNf, ierr, itime, itophys, itospec):
    """
    Convert data to be lmax & Nrhs
    """

    lvals = list(set(data[:,il])) # unique list of sorted lmax values
    lvals.sort()

    NRHS = 2*data[:,iNr]*data[:,iNf] # compute Nrhs

    uNRHS = list(set(NRHS)) # unique list of sorted Nrhs values
    uNRHS.sort()

    out = []

    for l in lvals:
        for Nrhs in uNRHS:

            # get all entries with this (\ell, Nrhs)
            ind = np.where((data[:,il]==l) & (NRHS==Nrhs))[0]

            if (len(ind) > 0): # found 1 or more valid entries
                entries = data[ind,:]

                avg = np.median(entries,axis=0) # representative value for all entries

                # store output
                entry = [l, Nrhs, avg[ierr], avg[itime], avg[itophys], avg[itospec]]
                out.append(entry)

    out = np.array(out) # (l,rhs,err,time,ToPhys,ToSpec)
    return out

def MatchData2(data, refdata, i1, i2):
    """
    Reorder entries in data to be same as refdata

    refdata = array of shape (Nruns, Ncols)
    i1,i2 are indices into Ncols axis that will be compared

    output is array same shape as refdata filled with entries from data
    """
    out = np.nan*np.zeros_like(refdata)

    Nruns = np.shape(refdata)[0]
    for i in range(Nruns):
        x = refdata[i,i1]
        y = refdata[i,i2]
        ind = np.where((x==data[:,i1]) & (y==data[:,i2]))
        if (len(ind[0]) == 1):
            out[i,:] = data[ind[0][0],:]
        elif (len(ind[0]) > 1):
            print("Multiple data entries: I1={}, I2={}".format(x,y))
    return out

def MatchData3(data, refdata, i1, i2, i3):
    """
    Reorder entries in data to be same as refdata

    refdata = array of shape (Nruns, Ncols)
    i1,i2,i3 are indices into Ncols axis that will be compared

    output is array same shape as refdata filled with entries from data
    """
    out = np.nan*np.zeros_like(refdata)

    Nruns = np.shape(refdata)[0]
    for i in range(Nruns):
        x = refdata[i,i1]
        y = refdata[i,i2]
        z = refdata[i,i3]
        ind = np.where((x==data[:,i1]) & \
                       (y==data[:,i2]) & \
                       (z==data[:,i3]))
        if (len(ind[0]) == 1):
            out[i,:] = data[ind[0][0],:]
        elif (len(ind[0]) > 1):
            print("Multiple data entries: I1={}, I2={}, I3={}".format(x,y,z))
    return out

def ConvertTo2D(x, y, z):
    # x,y,z = 1D arrays
    # output: X, Y, Z = 2D arrays

    # build unique 1D arrays for the grid
    xx = list(set(x))
    xx.sort()
    yy = list(set(y))
    yy.sort()

    # convert to integers
    xx = [int(i) for i in xx]
    yy = [int(i) for i in yy]

    # all will be shape (nx,ny)
    X, Y = np.meshgrid(xx, yy, indexing='ij')
    Z = np.zeros((len(xx), len(yy)))

    # fill in Z with proper data
    for i,_x in enumerate(xx):
        for j,_y in enumerate(yy):
            ind = np.where((x[:]==_x) & (y[:]==_y))
            if (len(ind[0]) == 1):
                Z[i,j] = z[ind[0][0]]
            else:
                Z[i,j] = np.nan
                if (len(ind[0]) == 0):
                    print("Missing data: x={}, y={}, setting to NaN".format(_x,_y))
                else:
                    print("Multiple entries: x={}, y={}, setting to NaN".format(_x,_y))

    return X, Y, Z

class MidPointNorm(Normalize):    
    """
    -----> Found this code on stackoverflow written by Annan <-----

    Center the colorbar by subclassing Normalize. To use it::

        norm = MidPointNorm(midpoint=3)
        imshow(X, norm=norm)

    """

    def __init__(self, midpoint=0, vmin=None, vmax=None, clip=False):
        Normalize.__init__(self, vmin, vmax, clip)
        self.midpoint = midpoint

    def __call__(self, value, clip=None):

        if clip is None:
            clip = self.clip

        result, is_scalar = self.process_value(value)

        self.autoscale_None(result)
        vmin, vmax, midpoint = self.vmin, self.vmax, self.midpoint

        if not (vmin < midpoint < vmax):
            raise ValueError("midpoint must be between maxvalue and minvalue.")       
        elif vmin == vmax:
            result.fill(0) # Or should it be all masked? Or 0.5?
        elif vmin > vmax:
            raise ValueError("maxvalue must be bigger than minvalue")
        else:
            vmin = float(vmin)
            vmax = float(vmax)
            if clip:
                mask = ma.getmask(result)
                result = ma.array(np.clip(result.filled(vmax), vmin, vmax),
                                  mask=mask)

            # ma division is very slow; we can take a shortcut
            resdat = result.data

            # First scale to -1 to 1 range, than to from 0 to 1.
            resdat -= midpoint            
            resdat[resdat>0] /= abs(vmax - midpoint)            
            resdat[resdat<0] /= abs(vmin - midpoint)

            resdat /= 2.
            resdat += 0.5
            result = ma.array(resdat, mask=result.mask, copy=False)                

        if is_scalar:
            result = result[0]            

        return result

    def inverse(self, value):

        if not self.scaled():
            raise ValueError("Not invertible until scaled")

        vmin, vmax, midpoint = self.vmin, self.vmax, self.midpoint

        if cbook.iterable(value):
            val = ma.asarray(value)
            val = 2 * (val-0.5)  
            val[val>0]  *= abs(vmax - midpoint)
            val[val<0] *= abs(vmin - midpoint)
            val += midpoint
            return val
        else:
            val = 2 * (val - 0.5)
            if val < 0: 
                return  val*abs(vmin-midpoint) + midpoint
            else:
                return  val*abs(vmax-midpoint) + midpoint

