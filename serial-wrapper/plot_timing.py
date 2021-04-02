"""
Plot Legendre transform timing between Rayleigh and SHTns

Usage:
    plot_timing.py [options] <SHT_timing_data> <Rayleigh_timing_data>

Options:
    --dpi=<d>        Resolution of images [default: 250]
    --output=<o>     Set the output basename for images [default: timing]
    --vs-lmax        Plot data versus l-max, default is N-rhs [default: False]
    --plot-all       Make individiual plots too [default: False]

"""

from __future__ import print_function
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.rcParams.update({"text.usetex":True})
import numpy as np

def compare(x,y): # simple comparison
    return np.allclose(x,y)

def convert_name(x): # convert Rayleigh & SHTns into integer codes
    x = str(x)
    if ("ray" in x.lower()):
       return 1
    elif ("sht" in x.lower()):
       return 0
    else:
       return -1

def MakeUnique(data, il, iNr, iNf, ierr, itime, itophys, itospec):

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

def main(Sdatafile, Rdatafile, dpi, basename, vs_lmax, all_plots):

    ilmax = 0; iNr = 1; iNf = 2; iNl = 3
    ierr = 4; itime = 5; itophys = 6; itospec = 7; icode = 8

    # (rows, columns)
    Rdata = np.loadtxt(Rdatafile, dtype=np.float64, comments="#",
                      delimiter=None, skiprows=0, usecols=None,
                      converters={icode:convert_name})
    # convert times to time/loop or time/call
    Rdata[:,itime]   = Rdata[:,itime]/Rdata[:,iNl]
    Rdata[:,itophys] = Rdata[:,itophys]/Rdata[:,iNl]
    Rdata[:,itospec] = Rdata[:,itospec]/Rdata[:,iNl]

    Sdata = np.loadtxt(Sdatafile, dtype=np.float64, comments="#",
                      delimiter=None, skiprows=0, usecols=None,
                      converters={icode:convert_name})
    # convert times to time/loop or time/call
    Sdata[:,itime]   = Sdata[:,itime]/Sdata[:,iNl]
    Sdata[:,itophys] = Sdata[:,itophys]/Sdata[:,iNl]
    Sdata[:,itospec] = Sdata[:,itospec]/Sdata[:,iNl]

    # convert variables from (lmax, Nr, Nf, err, times,...) to (lmax, Nrhs, err, times,...)
    Rdata = MakeUnique(Rdata, ilmax, iNr, iNf, ierr, itime, itophys, itospec)
    Sdata = MakeUnique(Sdata, ilmax, iNr, iNf, ierr, itime, itophys, itospec)

    # repoint indices
    ilmax = 0; iNrhs = 1; ierr = 2; itime = 3; itophys = 4; itospec = 5

    #######################################################

    # total time
    if (all_plots):
        ylabel = "Walltime (sec)"
        ys = [Rdata[:,itime], Sdata[:,itime]]
        labels = ["Rayleigh", "SHTns"]
        ls = ['-', ':']; cols = ['r', 'b']; marks = ['o', 'd']
        if (vs_lmax):
            xs = [Rdata[:,ilmax], Sdata[:,ilmax]]
            title  = r"Total Time vs $\ell_\mathrm{max}$"
            xlabel = r"$\ell_\mathrm{max}$"
            tag = "_total_time_vs_lmax.png"
        else:
            xs = [Rdata[:,iNrhs], Sdata[:,iNrhs]]
            title  = r"Total Time vs $N_\mathrm{rhs}$"
            xlabel = r"$N_\mathrm{rhs}$"
            tag = "_total_time_vs_Nrhs.png"

        output = basename + tag
        MakePlot(xs, ys, labels, title, xlabel, ylabel, output, ls, cols, marks, legend=True)

        # to-physical time
        ylabel = "Walltime (sec)"
        ys = [Rdata[:,itophys], Sdata[:,itophys]]
        labels = ["Rayleigh", "SHTns"]
        ls = ['-', ':']; cols = ['r', 'b']; marks = ['o', 'd']
        if (vs_lmax):
            xs = [Rdata[:,ilmax], Sdata[:,ilmax]]
            title  = r"To Physical Time vs $\ell_\mathrm{max}$"
            xlabel = r"$\ell_\mathrm{max}$"
            tag = "_tophysical_time_vs_lmax.png"
        else:
            xs = [Rdata[:,iNrhs], Sdata[:,iNrhs]]
            title  = r"To Physical Time vs $N_\mathrm{rhs}$"
            xlabel = r"$N_\mathrm{rhs}$"
            tag = "_tophysical_time_vs_Nrhs.png"

        output = basename + tag
        MakePlot(xs, ys, labels, title, xlabel, ylabel, output, ls, cols, marks, legend=True)

        # to-spectral time
        ylabel = "Walltime (sec)"
        ys = [Rdata[:,itospec], Sdata[:,itospec]]
        labels = ["Rayleigh", "SHTns"]
        ls = ['-', ':']; cols = ['r', 'b']; marks = ['o', 'd']
        if (vs_lmax):
            xs = [Rdata[:,ilmax], Sdata[:,ilmax]]
            title  = r"To Spectral Time vs $\ell_\mathrm{max}$"
            xlabel = r"$\ell_\mathrm{max}$"
            tag = "_tospectral_time_vs_lmax.png"
        else:
            xs = [Rdata[:,iNrhs], Sdata[:,iNrhs]]
            title  = r"To Spectral Time vs $N_\mathrm{rhs}$"
            xlabel = r"$N_\mathrm{rhs}$"
            tag = "_tospectral_time_vs_Nrhs.png"

        output = basename + tag
        MakePlot(xs, ys, labels, title, xlabel, ylabel, output, ls, cols, marks, legend=True)

    # errors
    ylabel = "Absolute Error"
    ys = [Rdata[:,ierr], Sdata[:,ierr]]
    labels = ["Rayleigh", "SHTns"]
    ls = ['-', ':']; cols = ['r', 'b']; marks = ['o', 'd']
    if (vs_lmax):
        xs = [Rdata[:,ilmax], Sdata[:,ilmax]]
        title  = r"Max Error vs $\ell_\mathrm{max}$"
        xlabel = r"$\ell_\mathrm{max}$"
        tag = "_error_vs_lmax.png"
    else:
        xs = [Rdata[:,iNrhs], Sdata[:,iNrhs]]
        title  = r"Max Error vs $N_\mathrm{rhs}$"
        xlabel = r"$N_\mathrm{rhs}$"
        tag = "_error_vs_Nrhs.png"

    output = basename + tag
    MakePlot(xs, ys, labels, title, xlabel, ylabel, output, ls, cols, marks, legend=True)

    # all times
    ylabel = "Walltime (sec)"
    ys = [Rdata[:,itime], Sdata[:,itime],
          Rdata[:,itophys], Sdata[:,itophys],
          Rdata[:,itospec], Sdata[:,itospec]]
    labels = ["Rayleigh Total", "SHTns Total",
              "Rayleigh ToPhys", "SHTns ToPhys",
              "Rayleigh ToSpec", "SHTns ToSpec"]
    ls = ['-', '-', ':', ':', '--', '--']; cols = ['r', 'b', 'r', 'b', 'r', 'b']
    marks = ['d', 'd', 'o', 'o', '*', '*']
    if (vs_lmax):
        xs = [Rdata[:,ilmax], Sdata[:,ilmax],
              Rdata[:,ilmax], Sdata[:,ilmax],
              Rdata[:,ilmax], Sdata[:,ilmax]]
        title  = r"Walltime vs $\ell_\mathrm{max}$"
        xlabel = r"$\ell_\mathrm{max}$"
        tag = "_times_vs_lmax.png"
    else:
        xs = [Rdata[:,iNrhs], Sdata[:,iNrhs],
              Rdata[:,iNrhs], Sdata[:,iNrhs],
              Rdata[:,iNrhs], Sdata[:,iNrhs]]
        title  = r"Walltime vs $N_\mathrm{rhs}$"
        xlabel = r"$N_\mathrm{rhs}$"
        tag = "_times_vs_Nrhs.png"

    output = basename + tag
    MakePlot(xs, ys, labels, title, xlabel, ylabel, output, ls, cols, marks, legend=True)

def MakePlot(xs, ys, labels, title, xlabel, ylabel, output,
             lstyles, colors, markers, legend=True, **plt_kw):

    fig = plt.figure(dpi=dpi)
    ax = fig.add_subplot(111)

    ax.set_title(title)

    #================================================================
    # axis setup
    tw = 1.5; th = 1.5
    font = 'large'; label = 'large'

    #----------------------------------------------------------------
    x = np.array(xs); xlim = (0.5*np.min(x), 2*np.max(x))
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
    y = np.array(ys); ylim = (0.5*np.min(y), 2*np.max(y))
    ax.set_ylim(ylim)
    ax.set_yscale('log')
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
    vs_lmax  = args['--vs-lmax']
    all_plots = args['--plot-all']

    main(Sdatafile, Rdatafile, dpi, output, vs_lmax, all_plots)

