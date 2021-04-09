from __future__ import print_function
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.rcParams.update({"text.usetex":True})
import numpy as np
from utils import ColorMap, MidPointNorm

def PlotGrid(x, y, data, output, xlabel, ylabel, title, xlog=True, ylog=True, dpi=250,
             **plt_kw):

    # for Pcolormesh, x & y need to be larger than data by 1 in each direction
    nx, ny = np.shape(data)
    X = np.zeros((nx+1,ny+1)) # add new row, new col, and new corner
    Y = np.zeros((nx+1,ny+1))
    xnew = 2*np.max(x); ynew = 2*np.max(y) # guess at "new" coordinates
    X[:-1,:-1] = x[:,:] # interior
    X[-1,:-1] = xnew    # last row
    X[:-1,-1] = x[:,-1] # last col
    X[-1,-1] = xnew     # new corner
    Y[:-1,:-1] = y[:,:] # interior
    Y[-1,:-1] = y[-1,:] # last row
    Y[:-1,-1] = ynew    # last col
    Y[-1,-1] = ynew     # new corner

    # mask the data, to get rid of NaNs
    data_masked = np.where(np.isfinite(data), data, -1.0)

    # choose colormap
    cmap = ColorMap('custom')
    cmap.set_under('k')       # values below vmin should be black

    fig = plt.figure(dpi=dpi)
    ax = fig.add_subplot(111)

    ax.set_title(title)
    
    #================================================================
    # axis setup
    tw = 1.5; th = 1.5
    font = 'large'; label = 'large'

    #----------------------------------------------------------------
    xlim = (0.5*np.min(X), 2*np.max(X))
    ax.set_xlim(xlim)
    if (xlog):
        ax.set_xscale('log')
    else:
        ax.set_xscale('linear')
    ax.set_xlabel(xlabel, fontsize=font)

    ax.tick_params(axis='x', which='both', direction='inout', labelsize=label,
                   width=tw, length=th)

    w = tw*mpl.rcParams['xtick.major.width']
    h = th*mpl.rcParams['xtick.major.size']
    ax.tick_params(axis='x', which='major', width=w, length=h)

    w = tw*mpl.rcParams['xtick.minor.width']
    h = th*mpl.rcParams['xtick.minor.size']
    ax.tick_params(axis='x', which='minor', width=w, length=h)

    #----------------------------------------------------------------
    ylim = (0.5*np.min(Y), 2*np.max(Y))
    ax.set_ylim(ylim)
    if (ylog):
        ax.set_yscale('log')
    else:
        ax.set_yscale('linear')
    ax.set_ylabel(ylabel, fontsize=font)

    ax.tick_params(axis='y', which='both', direction='inout', labelsize=label,
                   width=tw, length=th)

    w = tw*mpl.rcParams['ytick.major.width']
    h = th*mpl.rcParams['ytick.major.size']
    ax.tick_params(axis='x', which='major', width=w, length=h)

    w = tw*mpl.rcParams['ytick.minor.width']
    h = th*mpl.rcParams['ytick.minor.size']
    ax.tick_params(axis='y', which='minor', width=w, length=h)

    #================================================================
    vmin = np.min(data_masked[np.where(data_masked>0.0)]); vmax = np.max(data_masked)
    if ((vmin < 1.0) and (1.0 < vmax)):
        mid = 1.0
    else:
        mid = 0.5*(vmax+vmin)
    img = plt.pcolormesh(X, Y, data_masked, cmap=cmap, vmin=vmin, vmax=vmax,
                         shading='nearest', norm=MidPointNorm(midpoint=mid))

    # colorbar
    cbar = fig.colorbar(img, orientation='vertical')
    cbar.ax.tick_params(direction='out', labelsize=label)

    #================================================================
    plt.tight_layout()

    plt.savefig(output, bbox_inches='tight', dpi=dpi)
    print("\tsaved image: {}\n".format(output))

    plt.clf(); plt.close() # clear/close figure when done

    return

def MakePlot(xs, ys, labels, title, xlabel, ylabel, output,
             lstyles, colors, markers, legend=True, ylog=True, dpi=250, **plt_kw):

    fig = plt.figure(dpi=dpi)
    ax = fig.add_subplot(111)

    ax.set_title(title)

    #================================================================
    # axis setup
    tw = 1.5; th = 1.5
    font = 'large'; label = 'large'

    #----------------------------------------------------------------
    xmin = np.min([np.min(i) for i in xs])
    xmax = np.max([np.max(i) for i in xs])
    xlim = (0.5*xmin, 2*xmax)
    ax.set_xlim(xlim)
    ax.set_xscale('log')
    ax.set_xlabel(xlabel, fontsize=font)

    ax.tick_params(axis='x', which='both', direction='inout', labelsize=label,
                   width=tw, length=th)

    w = tw*mpl.rcParams['xtick.major.width']
    h = th*mpl.rcParams['xtick.major.size']
    ax.tick_params(axis='x', which='major', width=w, length=h)

    w = tw*mpl.rcParams['xtick.minor.width']
    h = th*mpl.rcParams['xtick.minor.size']
    ax.tick_params(axis='x', which='minor', width=w, length=h)

    #----------------------------------------------------------------
    ymin = np.min([np.min(i) for i in ys])
    ymax = np.max([np.max(i) for i in ys])
    ylim = (0.5*ymin, 2*ymax)
    ax.set_ylim(ylim)
    if (ylog):
        ax.set_yscale('log')
    else:
        ax.set_yscale('linear')
    ax.set_ylabel(ylabel, fontsize=font)

    ax.tick_params(axis='y', which='both', direction='inout', labelsize=label,
                   width=tw, length=th)

    w = tw*mpl.rcParams['ytick.major.width']
    h = th*mpl.rcParams['ytick.major.size']
    ax.tick_params(axis='x', which='major', width=w, length=h)

    w = tw*mpl.rcParams['ytick.minor.width']
    h = th*mpl.rcParams['ytick.minor.size']
    ax.tick_params(axis='y', which='minor', width=w, length=h)

    #================================================================
    for x,y,l,ls,c,m in zip(xs, ys, labels, lstyles, colors, markers):
        ax.plot(x, y, label=l, marker=m, ls=ls, color=c, **plt_kw)

    #================================================================
    if (legend):
        plt.legend(loc='best', fontsize='large', ncol=1, numpoints=1,
               frameon=True, framealpha=0.5)

    #================================================================
    plt.tight_layout()

    plt.savefig(output, bbox_inches='tight', dpi=dpi)
    print("\tsaved image: {}\n".format(output))

    plt.clf(); plt.close() # clear/close figure when done

