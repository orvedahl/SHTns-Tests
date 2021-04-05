"""
Plot Legendre transform timing between Rayleigh and SHTns

Usage:
    plot_timing.py [options] <Rayleigh_timing_data> <SHT_timing_data>

Options:
    --dpi=<d>     Resolution of images [default: 250]
    --output=<o>  Set the output basename for images [default: timing]
    --lmax=<l>    Plot data vs Nrhs for the given lmax value
    --Nrhs=<n>    Plot data vs lmax for the given Nrhs value
    --plot-all    Make individiual plots too [default: False]
    --grid        Plot speedup data in a grid [default: False]

"""

from __future__ import print_function
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.rcParams.update({"text.usetex":True})
import numpy as np
from utils import MatchData2, MakeUnique, MidPointNorm, ConvertTo2D, ColorMap

def convert_name(x): # convert Rayleigh & SHTns into integer codes
    x = str(x)
    if ("ray" in x.lower()):
       return 1
    elif ("sht" in x.lower()):
       return 0
    else:
       return -1

def main(Sdatafile, Rdatafile, dpi, basename, lmax_choose, Nrhs_choose, all_plots, grid):

    if ((lmax_choose is None) and (Nrhs_choose is None) and (not grid)):
        print("\nPick something to do: either --lmax, --Nrhs, or --grid\n")
        return

    ilmax = 0; iNr = 1; iNf = 2; iNl = 3
    ierr = 4; itime = 5; itophys = 6; itospec = 7; icode = 8

    # (rows, columns)
    Rdata_o = np.loadtxt(Rdatafile, dtype=np.float64, comments="#",
                      delimiter=None, skiprows=0, usecols=None,
                      converters={icode:convert_name})
    # convert times to time/loop or time/call
    Rdata_o[:,itime]   = Rdata_o[:,itime]/Rdata_o[:,iNl]
    Rdata_o[:,itophys] = Rdata_o[:,itophys]/Rdata_o[:,iNl]
    Rdata_o[:,itospec] = Rdata_o[:,itospec]/Rdata_o[:,iNl]
    Rdata_o[:,ierr]    = Rdata_o[:,ierr]/Rdata_o[:,iNl]

    Sdata_o = np.loadtxt(Sdatafile, dtype=np.float64, comments="#",
                      delimiter=None, skiprows=0, usecols=None,
                      converters={icode:convert_name})
    # convert times to time/loop or time/call
    Sdata_o[:,itime]   = Sdata_o[:,itime]/Sdata_o[:,iNl]
    Sdata_o[:,itophys] = Sdata_o[:,itophys]/Sdata_o[:,iNl]
    Sdata_o[:,itospec] = Sdata_o[:,itospec]/Sdata_o[:,iNl]
    Sdata_o[:,ierr]    = Sdata_o[:,ierr]/Sdata_o[:,iNl]

    # convert from (lmax, Nr, Nf, err, times,...) to (lmax, Nrhs, err, times,...)
    Rdata = MakeUnique(Rdata_o, ilmax, iNr, iNf, ierr, itime, itophys, itospec)
    Sdata = MakeUnique(Sdata_o, ilmax, iNr, iNf, ierr, itime, itophys, itospec)

    # repoint indices
    ilmax = 0; iNrhs = 1; ierr = 2; itime = 3; itophys = 4; itospec = 5

    # make Sdata & Rdata the same dimensions and ordering (NaN for mismatched points)
    n1 = np.shape(Sdata)[0]; n2 = np.shape(Rdata)[0]
    if (n1 > n2):
        # Sdata has more runs, use it as reference
        refdata = Sdata
        Rdata = MatchData2(Rdata, refdata, ilmax, iNrhs)
    else:
        # Rdata has more runs, use it as reference
        refdata = Rdata
        Sdata = MatchData2(Sdata, refdata, ilmax, iNrhs)

    lvals = list(set(refdata[:,ilmax])); lvals.sort()
    Nrhs = list(set(refdata[:,iNrhs])); Nrhs.sort()
    lvals = [int(i) for i in lvals]
    Nrhs = [int(i) for i in Nrhs]

    #######################################################

    if (grid): # plot grid of results
        x = refdata[:,ilmax]; xlabel = r"$\ell_\mathrm{max}$"
        y = refdata[:,iNrhs]; ylabel = r"$N_\mathrm{rhs}$"

        X, Y, Z = ConvertTo2D(x, y, Rdata[:,itime]/Sdata[:,itime])
        output = basename + "_total_speedup_grid.png"
        PlotGrid(X, Y, Z, output, xlabel, ylabel, "Total Time: SHTns Speedup Compared to Rayleigh")

        X, Y, Z = ConvertTo2D(x, y, Rdata[:,itophys]/Sdata[:,itophys])
        output = basename + "_tophys_speedup_grid.png"
        PlotGrid(X, Y, Z, output, xlabel, ylabel, "ToPhys: SHTns Speedup Compared to Rayleigh")

        X, Y, Z = ConvertTo2D(x, y, Rdata[:,itospec]/Sdata[:,itospec])
        output = basename + "_tospec_speedup_grid.png"
        PlotGrid(X, Y, Z, output, xlabel, ylabel, "ToSpec: SHTns Speedup Compared to Rayleigh")

    if (lmax_choose is not None): # plot data vs Nrhs for the chosen l_max
        lmax_choose = int(lmax_choose)
        if (lmax_choose not in lvals):
            print("\nChosen value of l_max does not appear in the data: lmax={}".format(lmax_choose))
        else:

            ind = np.where(refdata[:,ilmax]==lmax_choose)[0]
            Rd = Rdata[ind,:]; Sd = Sdata[ind,:]

            xlabel = r"$N_\mathrm{rhs}$"

            if (all_plots):
                # total time
                title  = r"Total Time vs $N_\mathrm{rhs}$"
                ylabel = "Walltime/loop (sec)"
                xs = [Rd[:,iNrhs], Sd[:,iNrhs]]
                ys = [Rd[:,itime], Sd[:,itime]]
                labels = ["Rayleigh", "SHTns"]
                ls = ['-', ':']; cols = ['r', 'b']; marks = ['o', 'd']
                output = basename + "_total_time_vs_Nrhs_lmax{}.png".format(lmax_choose)
                MakePlot(xs, ys, labels, title, xlabel, ylabel, output, ls, cols, marks, legend=True)

                # to-physical time
                title  = r"To Physical Time vs $N_\mathrm{rhs}$"
                ylabel = "Walltime/loop (sec)"
                xs = [Rd[:,iNrhs], Sd[:,iNrhs]]
                ys = [Rd[:,itophys], Sd[:,itophys]]
                labels = ["Rayleigh", "SHTns"]
                ls = ['-', ':']; cols = ['r', 'b']; marks = ['o', 'd']
                output = basename + "_tophysical_time_vs_Nrhs_lmax{}.png".format(lmax_choose)
                MakePlot(xs, ys, labels, title, xlabel, ylabel, output, ls, cols, marks, legend=True)

                # to-spectral time
                title  = r"To Spectral Time vs $N_\mathrm{rhs}$"
                ylabel = "Walltime/loop (sec)"
                xs = [Rd[:,iNrhs], Sd[:,iNrhs]]
                ys = [Rd[:,itospec], Sd[:,itospec]]
                labels = ["Rayleigh", "SHTns"]
                ls = ['-', ':']; cols = ['r', 'b']; marks = ['o', 'd']
                output = basename + "_tospectral_time_vs_Nrhs_lmax{}.png".format(lmax_choose)
                MakePlot(xs, ys, labels, title, xlabel, ylabel, output, ls, cols, marks, legend=True)

            # errors
            title  = r"Max Error vs $N_\mathrm{rhs}$"
            ylabel = "Absolute Error/loop"
            xs = [Rd[:,iNrhs], Sd[:,iNrhs]]
            ys = [Rd[:,ierr], Sd[:,ierr]]
            labels = ["Rayleigh", "SHTns"]
            ls = ['-', ':']; cols = ['r', 'b']; marks = ['o', 'd']
            output = basename + "_error_vs_Nrhs_lmax{}.png".format(lmax_choose)
            MakePlot(xs, ys, labels, title, xlabel, ylabel, output, ls, cols, marks, legend=True)

            # all times
            title  = r"Walltime vs $N_\mathrm{rhs}$"
            ylabel = "Walltime/loop (sec)"
            xs = [Rd[:,iNrhs], Sd[:,iNrhs],
                  Rd[:,iNrhs], Sd[:,iNrhs],
                  Rd[:,iNrhs], Sd[:,iNrhs]]
            ys = [Rd[:,itime], Sd[:,itime],
                  Rd[:,itophys], Sd[:,itophys],
                  Rd[:,itospec], Sd[:,itospec]]
            labels = ["Rayleigh Total", "SHTns Total",
                      "Rayleigh ToPhys", "SHTns ToPhys",
                      "Rayleigh ToSpec", "SHTns ToSpec"]
            ls = ['-', '-', ':', ':', '--', '--']; cols = ['r', 'b', 'r', 'b', 'r', 'b']
            marks = ['d', 'd', 'o', 'o', '*', '*']
            xlabel = r"$N_\mathrm{rhs}$"
            output = basename + "_times_vs_Nrhs_lmax{}.png".format(lmax_choose)
            MakePlot(xs, ys, labels, title, xlabel, ylabel, output, ls, cols, marks, legend=True)

    if (Nrhs_choose is not None): # plot data vs lmax for the chosen Nrhs
        Nrhs_choose = int(Nrhs_choose)
        if (Nrhs_choose not in Nrhs):
            print("\nChosen value of Nrhs does not appear in the data: Nrhs={}".format(Nrhs_choose))
        else:

            ind = np.where(refdata[:,iNrhs]==Nrhs_choose)[0]
            Rd = Rdata[ind,:]; Sd = Sdata[ind,:]

            xlabel = r"$\ell_\mathrm{max}$"

            if (all_plots):
                # total time
                title  = r"Total Time vs $\ell_\mathrm{max}$"
                ylabel = "Walltime/loop (sec)"
                xs = [Rd[:,ilmax], Sd[:,ilmax]]
                ys = [Rd[:,itime], Sd[:,itime]]
                labels = ["Rayleigh", "SHTns"]
                ls = ['-', ':']; cols = ['r', 'b']; marks = ['o', 'd']
                output = basename + "_total_time_vs_lmax_Nrhs{}.png".format(Nrhs_choose)
                MakePlot(xs, ys, labels, title, xlabel, ylabel, output, ls, cols, marks, legend=True)

                # to-physical time
                title  = r"To Physical Time vs $\ell_\mathrm{max}$"
                ylabel = "Walltime/loop (sec)"
                xs = [Rd[:,ilmax], Sd[:,ilmax]]
                ys = [Rd[:,itophys], Sd[:,itophys]]
                labels = ["Rayleigh", "SHTns"]
                ls = ['-', ':']; cols = ['r', 'b']; marks = ['o', 'd']
                output = basename + "_tophysical_time_vs_lmax_Nrhs{}.png".format(Nrhs_choose)
                MakePlot(xs, ys, labels, title, xlabel, ylabel, output, ls, cols, marks, legend=True)

                # to-spectral time
                title  = r"To Spectral Time vs $\ell_\mathrm{max}$"
                ylabel = "Walltime/loop (sec)"
                xs = [Rd[:,ilmax], Sd[:,ilmax]]
                ys = [Rd[:,itospec], Sd[:,itospec]]
                labels = ["Rayleigh", "SHTns"]
                ls = ['-', ':']; cols = ['r', 'b']; marks = ['o', 'd']
                output = basename + "_tospectral_time_vs_lmax_Nrhs{}.png".format(Nrhs_choose)
                MakePlot(xs, ys, labels, title, xlabel, ylabel, output, ls, cols, marks, legend=True)

            # errors
            title  = r"Max Error vs $\ell_\mathrm{max}$"
            ylabel = "Absolute Error/loop"
            xs = [Rd[:,ilmax], Sd[:,ilmax]]
            ys = [Rd[:,ierr], Sd[:,ierr]]
            labels = ["Rayleigh", "SHTns"]
            ls = ['-', ':']; cols = ['r', 'b']; marks = ['o', 'd']
            output = basename + "_error_vs_lmax_Nrhs{}.png".format(Nrhs_choose)
            MakePlot(xs, ys, labels, title, xlabel, ylabel, output, ls, cols, marks, legend=True)

            # all times
            title  = r"Walltime vs $\ell_\mathrm{max}$"
            ylabel = "Walltime/loop (sec)"
            xs = [Rd[:,ilmax], Sd[:,ilmax],
                  Rd[:,ilmax], Sd[:,ilmax],
                  Rd[:,ilmax], Sd[:,ilmax]]
            ys = [Rd[:,itime], Sd[:,itime],
                  Rd[:,itophys], Sd[:,itophys],
                  Rd[:,itospec], Sd[:,itospec]]
            labels = ["Rayleigh Total", "SHTns Total",
                      "Rayleigh ToPhys", "SHTns ToPhys",
                      "Rayleigh ToSpec", "SHTns ToSpec"]
            ls = ['-', '-', ':', ':', '--', '--']; cols = ['r', 'b', 'r', 'b', 'r', 'b']
            marks = ['d', 'd', 'o', 'o', '*', '*']
            output = basename + "_times_vs_lmax_Nrhs{}.png".format(Nrhs_choose)
            MakePlot(xs, ys, labels, title, xlabel, ylabel, output, ls, cols, marks, legend=True)

def PlotGrid(x, y, data, output, xlabel, ylabel, title, xlog=True, ylog=True, **plt_kw):

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
             lstyles, colors, markers, legend=True, ylog=True, **plt_kw):

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

if __name__ == "__main__":

    from docopt import docopt
    args = docopt(__doc__)

    Sdatafile = args['<SHT_timing_data>']
    Rdatafile = args['<Rayleigh_timing_data>']
    dpi      = float(args['--dpi'])
    output   = args['--output']
    lmax_choose = args['--lmax']
    Nrhs_choose = args['--Nrhs']
    all_plots = args['--plot-all']
    grid = args['--grid']

    main(Sdatafile, Rdatafile, dpi, output, lmax_choose, Nrhs_choose, all_plots, grid)

