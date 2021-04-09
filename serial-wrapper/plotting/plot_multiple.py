"""
Plot Legendre transform timing between Rayleigh and SHTns

Usage:
    plot_timing.py [options]

Options:
    --dpi=<d>       Resolution of images [default: 250]
    --output=<o>    Set the output basename for images [default: timing]
    --lmax=<l>      Plot data vs Nrhs for the given lmax value
    --Nrhs=<n>      Plot data vs lmax for the given Nrhs value
    --plot-all      Make individiual plots too [default: False]
    --grid          Plot speedup data in a grid [default: False]
    --rayleigh=<r>  List of Rayleigh files, accepts globbing characters
    --shtns=<s>     List of SHTns files, accepts globbing characters
"""

from __future__ import print_function
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.rcParams.update({"text.usetex":True})
import numpy as np
from utils import MatchData2, MakeUnique, ConvertTo2D
from plotting import PlotGrid, MakePlot

def convert_name(x): # convert Rayleigh & SHTns into integer codes
    x = str(x)
    if ("ray" in x.lower()):
       return 1
    elif ("sht" in x.lower()):
       return 0
    else:
       return -1

def main(Sfiles, Rfiles, dpi, basename, lmax_choose, Nrhs_choose, all_plots, grid):

    if ((lmax_choose is None) and (Nrhs_choose is None) and (not grid)):
        print("\nPick something to do: either --lmax, --Nrhs, or --grid\n")
        return

    ilmax = 0; iNr = 1; iNf = 2; iNl = 3
    ierr = 4; itime = 5; itophys = 6; itospec = 7; icode = 8

    if ((len(Sfiles) < 1) or (len(Rfiles) < 1)):
        if (len(Sfiles) < 1):
            print("ERROR: must provide at least one SHTns datafile")
        if (len(Rfiles) < 1):
            print("ERROR: must provide at least one Rayleigh datafile")
        return

    load_args = dict(dtype=np.float64, comments="#",
                     delimiter=None, skiprows=0, usecols=None,
                     converters={icode:convert_name})

    Rdata = {}
    for n in range(len(Rfiles)):
        # (rows, columns)
        Rdata[n] = np.loadtxt(Rfiles[n], **load_args)

        # convert times to time/loop or time/call
        Rdata[n][:,itime]   = Rdata[n][:,itime]/Rdata[n][:,iNl]
        Rdata[n][:,itophys] = Rdata[n][:,itophys]/Rdata[n][:,iNl]
        Rdata[n][:,itospec] = Rdata[n][:,itospec]/Rdata[n][:,iNl]
        Rdata[n][:,ierr]    = Rdata[n][:,ierr]/Rdata[n][:,iNl]

    Sdata = {}
    for n in range(len(Sfiles)):
        # (rows, columns)
        Sdata[n] = np.loadtxt(Sfiles[n], **load_args)

        # convert times to time/loop or time/call
        Sdata[n][:,itime]   = Sdata[n][:,itime]/Sdata[n][:,iNl]
        Sdata[n][:,itophys] = Sdata[n][:,itophys]/Sdata[n][:,iNl]
        Sdata[n][:,itospec] = Sdata[n][:,itospec]/Sdata[n][:,iNl]
        Sdata[n][:,ierr]    = Sdata[n][:,ierr]/Sdata[n][:,iNl]

    # convert from (lmax, Nr, Nf, err, times,...) to (lmax, Nrhs, err, times,...)
    nRruns = 0
    for n in Rdata.keys():
        Rdata[n] = MakeUnique(Rdata[n], ilmax, iNr, iNf, ierr, itime, itophys, itospec)

        N = np.shape(Rdata[n])[0]
        if (N > nRruns): # find key to dataset with largest number of runs
            nRruns = N
            iRkey = n

    nSruns = 0
    for n in Sdata.keys():
        Sdata[n] = MakeUnique(Sdata[n], ilmax, iNr, iNf, ierr, itime, itophys, itospec)

        N = np.shape(Sdata[n])[0]
        if (N > nSruns): # find key to dataset with largest number of runs
            nSruns = N
            iSkey = n

    # repoint indices
    ilmax = 0; iNrhs = 1; ierr = 2; itime = 3; itophys = 4; itospec = 5

    # make all Sdata datasets the same size/order
    Srefdata = Sdata[iSkey]
    for n in Sdata.keys():
        if (n == iSkey): continue
        Sdata[n] = MatchData2(Sdata[n], Srefdata, ilmax, iNrhs)

    # make all Rdata datasets the same size/order
    Rrefdata = Rdata[iRkey]
    for n in Rdata.keys():
        if (n == iRkey): continue
        Rdata[n] = MatchData2(Rdata[n], Rrefdata, ilmax, iNrhs)

    # convert dictionaries into arrays
    shp = np.shape(Sdata[iSkey])
    Sdatasets = np.zeros((len(Sdata.keys()),) + shp)
    for i,n in enumerate(Sdata.keys()):
        Sdatasets[i,...] = Sdata[n][...]

    shp = np.shape(Rdata[iRkey])
    Rdatasets = np.zeros((len(Rdata.keys()),) + shp)
    for i,n in enumerate(Rdata.keys()):
        Rdatasets[i,...] = Rdata[n][...]

    # compute averages over all datasets, use masked arrays to account for possible NaNs
    Sexclude = np.logical_not(np.isfinite(Sdatasets))
    Sdata = np.ma.array(Sdatasets, mask=Sexclude)
    Sdata = Sdata.mean(axis=0)

    Rexclude = np.logical_not(np.isfinite(Rdatasets))
    Rdata = np.ma.array(Rdatasets, mask=Rexclude)
    Rdata = Rdata.mean(axis=0)

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
        PlotGrid(X, Y, Z, output, xlabel, ylabel, "Total Time: SHTns Speedup Compared to Rayleigh", dpi=dpi)

        X, Y, Z = ConvertTo2D(x, y, Rdata[:,itophys]/Sdata[:,itophys])
        output = basename + "_tophys_speedup_grid.png"
        PlotGrid(X, Y, Z, output, xlabel, ylabel, "ToPhys: SHTns Speedup Compared to Rayleigh", dpi=dpi)

        X, Y, Z = ConvertTo2D(x, y, Rdata[:,itospec]/Sdata[:,itospec])
        output = basename + "_tospec_speedup_grid.png"
        PlotGrid(X, Y, Z, output, xlabel, ylabel, "ToSpec: SHTns Speedup Compared to Rayleigh", dpi=dpi)

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
                MakePlot(xs, ys, labels, title, xlabel, ylabel, output, ls, cols, marks, legend=True, dpi=dpi)

                # to-physical time
                title  = r"To Physical Time vs $N_\mathrm{rhs}$"
                ylabel = "Walltime/loop (sec)"
                xs = [Rd[:,iNrhs], Sd[:,iNrhs]]
                ys = [Rd[:,itophys], Sd[:,itophys]]
                labels = ["Rayleigh", "SHTns"]
                ls = ['-', ':']; cols = ['r', 'b']; marks = ['o', 'd']
                output = basename + "_tophysical_time_vs_Nrhs_lmax{}.png".format(lmax_choose)
                MakePlot(xs, ys, labels, title, xlabel, ylabel, output, ls, cols, marks, legend=True, dpi=dpi)

                # to-spectral time
                title  = r"To Spectral Time vs $N_\mathrm{rhs}$"
                ylabel = "Walltime/loop (sec)"
                xs = [Rd[:,iNrhs], Sd[:,iNrhs]]
                ys = [Rd[:,itospec], Sd[:,itospec]]
                labels = ["Rayleigh", "SHTns"]
                ls = ['-', ':']; cols = ['r', 'b']; marks = ['o', 'd']
                output = basename + "_tospectral_time_vs_Nrhs_lmax{}.png".format(lmax_choose)
                MakePlot(xs, ys, labels, title, xlabel, ylabel, output, ls, cols, marks, legend=True, dpi=dpi)

            # errors
            title  = r"Max Error vs $N_\mathrm{rhs}$"
            ylabel = "Absolute Error/loop"
            xs = [Rd[:,iNrhs], Sd[:,iNrhs]]
            ys = [Rd[:,ierr], Sd[:,ierr]]
            labels = ["Rayleigh", "SHTns"]
            ls = ['-', ':']; cols = ['r', 'b']; marks = ['o', 'd']
            output = basename + "_error_vs_Nrhs_lmax{}.png".format(lmax_choose)
            MakePlot(xs, ys, labels, title, xlabel, ylabel, output, ls, cols, marks, legend=True, dpi=dpi)

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
            MakePlot(xs, ys, labels, title, xlabel, ylabel, output, ls, cols, marks, legend=True, dpi=dpi)

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
                MakePlot(xs, ys, labels, title, xlabel, ylabel, output, ls, cols, marks, legend=True, dpi=dpi)

                # to-physical time
                title  = r"To Physical Time vs $\ell_\mathrm{max}$"
                ylabel = "Walltime/loop (sec)"
                xs = [Rd[:,ilmax], Sd[:,ilmax]]
                ys = [Rd[:,itophys], Sd[:,itophys]]
                labels = ["Rayleigh", "SHTns"]
                ls = ['-', ':']; cols = ['r', 'b']; marks = ['o', 'd']
                output = basename + "_tophysical_time_vs_lmax_Nrhs{}.png".format(Nrhs_choose)
                MakePlot(xs, ys, labels, title, xlabel, ylabel, output, ls, cols, marks, legend=True, dpi=dpi)

                # to-spectral time
                title  = r"To Spectral Time vs $\ell_\mathrm{max}$"
                ylabel = "Walltime/loop (sec)"
                xs = [Rd[:,ilmax], Sd[:,ilmax]]
                ys = [Rd[:,itospec], Sd[:,itospec]]
                labels = ["Rayleigh", "SHTns"]
                ls = ['-', ':']; cols = ['r', 'b']; marks = ['o', 'd']
                output = basename + "_tospectral_time_vs_lmax_Nrhs{}.png".format(Nrhs_choose)
                MakePlot(xs, ys, labels, title, xlabel, ylabel, output, ls, cols, marks, legend=True, dpi=dpi)

            # errors
            title  = r"Max Error vs $\ell_\mathrm{max}$"
            ylabel = "Absolute Error/loop"
            xs = [Rd[:,ilmax], Sd[:,ilmax]]
            ys = [Rd[:,ierr], Sd[:,ierr]]
            labels = ["Rayleigh", "SHTns"]
            ls = ['-', ':']; cols = ['r', 'b']; marks = ['o', 'd']
            output = basename + "_error_vs_lmax_Nrhs{}.png".format(Nrhs_choose)
            MakePlot(xs, ys, labels, title, xlabel, ylabel, output, ls, cols, marks, legend=True, dpi=dpi)

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
            MakePlot(xs, ys, labels, title, xlabel, ylabel, output, ls, cols, marks, legend=True, dpi=dpi)

if __name__ == "__main__":

    from docopt import docopt
    from glob import glob

    args = docopt(__doc__)

    Sfiles = args['--shtns']
    Rfiles = args['--rayleigh']
    dpi      = float(args['--dpi'])
    output   = args['--output']
    lmax_choose = args['--lmax']
    Nrhs_choose = args['--Nrhs']
    all_plots = args['--plot-all']
    grid = args['--grid']

    # convert --shtns & --rayleigh to list of files
    if (Sfiles is not None):
       Sfiles = glob(Sfiles)
    else:
       raise ValueError("Must provide SHTns files with --shtns=...")

    if (Rfiles is not None):
       Rfiles = glob(Rfiles)
    else:
       raise ValueError("Must provide Rayleigh files with --rayleigh=...")

    main(Sfiles, Rfiles, dpi, output, lmax_choose, Nrhs_choose, all_plots, grid)

