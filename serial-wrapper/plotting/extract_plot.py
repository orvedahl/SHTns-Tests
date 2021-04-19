"""
Extract performance data from a collection of files and plot it

Usage:
    extract.py [options]

Options:
    --r-files=<r>  Collection of pure Rayleigh files to process
    --s-files=<s>  Collection of Rayleigh-using-SHTns files to process
    --m-files=<m>  Collection of MagIC files to process
    --output=<o>   Set the output filename [default: scaling.png]
    --dpi=<d>      Resolution of image [default: 250]
"""

from __future__ import print_function
import matplotlib.pyplot as plt
plt.rcParams.update({"text.usetex":True})
import numpy as np
from plotting import MakePlot
from utils import ColorMarks

def ParseMagic(f):
    """
    parse the standard out of a MagIC simulation for scaling data
    returns: lmax, Nr, iter/sec, elapsed time, Ncpu, Nprow, Npcol
    """
    raise ValueError("MagIC input files not implemented yet")
    return lmax, nr, iterrate, time, ncpu, nprow, npcol

def ParseRayleigh(f):
    """
    parse the standard out of a Rayleigh simulation for scaling data
    returns: lmax, Nr, iter/sec, elapsed time, Ncpu, Nprow, Npcol
    """
    lmax = nr = iterrate = time = ncpu = nprow = npcol = None
    with open(f, "r") as mf:
        for line in mf:
            if ("ncpu" in line.lower()):
                l = line.split(":")[1]
                ncpu = int(l.strip())
            elif ("nprow" in line.lower()):
                l = line.split(":")[1]
                nprow = int(l.strip())
            elif ("npcol" in line.lower()):
                l = line.split(":")[1]
                npcol = int(l.strip())
            elif ("n_r" in line.lower()):
                l = line.split(":")[1]
                nr = int(l.strip())
            elif ("ell_max" in line.lower()):
                l = line.split(":")[1]
                lmax = int(l.strip())
            elif ("elapsed" in line.lower()):
                l = line.split(":")[1]
                time = float(l.strip())
            elif (("iteration" not in line.lower()) and ("iter/sec" in line.lower())):
                l = line.split(":")[1]
                iterrate = float(l.strip())

    return lmax, nr, iterrate, time, ncpu, nprow, npcol

def ParseFiles(files, magic=False):

    if (files is None): return None

    if (magic): # define interface
        Parse = ParseMagic
    else:
        Parse = ParseRayleigh

    results = []
    for f in files:
        res = Parse(f)
        if (None in res): continue # case did not succeed or did not finish
        results.append(res)

    return np.array(results)

def SplitIntoRuns(data, inds):
    """
    Reorganize data into collection of (run vs core) entries
    data is (N_all_runs, Nentries), inds is dictionary of indices.

    returns:
        -list of (Nrun, Nentries) arrays at a particular [lmax, nr]
        -list of [lmax, Nr] items for each array
    """
    if (data is None): return [], []

    ir = inds['nr']; il = inds['lmax']; ic = inds['cores']
    all_nr = data[:,ir]
    all_lmax = data[:,il]

    # make nr/lmax unique and sorted
    all_nr = list(set(all_nr))
    all_lmax = list(set(all_lmax))

    all_nr.sort()
    all_lmax.sort()

    results = []
    resolution = []
    for l in all_lmax: # loop over unique lmax/nr
        for r in all_nr:
            I = np.where((data[:,ir]==r) & (data[:,il]==l))[0]

            if (len(I) > 0): # found one or more entries at this (lmax,nr)

                entry = data[I,:] # get all suitable entries

                i = np.argsort(entry[:,ic]) # sort entry based on core count and store
                results.append(entry[i])
                resolution.append([l, r])

    return results, resolution

def main(Rfiles, Sfiles, Mfiles, output, dpi):

    test = [x is None for x in [Rfiles,Sfiles,Mfiles]]
    if (all(test)):
        print("\nMust provide at least one of --r-files, --s-files, --m-files\n")
        return

    # set data and indices
    ind = {'lmax':0, 'nr':1, 'iter':2, 'time':3, 'cores':4, 'nprow':5, 'npcol':6}
    Rdata = ParseFiles(Rfiles) # (Nruns, Nentry)
    Sdata = ParseFiles(Sfiles) # (Nruns, Nentry)
    Mdata = ParseFiles(Mfiles, magic=True) # (Nruns, Nentry)

    # split into unique runs, i.e., same resolution but different core count
    Rruns, Rres = SplitIntoRuns(Rdata, ind)
    Sruns, Sres = SplitIntoRuns(Sdata, ind)
    Mruns, Mres = SplitIntoRuns(Mdata, ind)

    # make plot
    xlabel = r"$N_\mathrm{cores}$"
    ylabel = "Walltime (sec)"
    title = "Scaling Results"

    xs = []; ys = []; labels = []; ls = []; colors = []; markers = []

    color_marks = ColorMarks()
    ic = ind['cores']; it = ind['time']
    for i in range(len(Rruns)):
        x = Rruns[i][:,ic]
        y = Rruns[i][:,it]

        lmax = int(Rres[i][0]); nr = int(Rres[i][1])
        l = r"Rayleigh, $\ell_\mathrm{{max}}={{{}}}$ $N_r={{{}}}$".format(lmax,nr)

        lstyle = '-'
        c, m = color_marks('slow-color')

        xs.append(x); ys.append(y)
        ls.append(lstyle); colors.append(c); markers.append(m); labels.append(l)

    color_marks.reset_counters(); color_marks.cind += 1
    for i in range(len(Sruns)):
        x = Sruns[i][:,ic]
        y = Sruns[i][:,it]

        lmax = int(Sres[i][0]); nr = int(Sres[i][1])
        l = r"SHTns, $\ell_\mathrm{{max}}={{{}}}$ $N_r={{{}}}$".format(lmax,nr)

        lstyle = '-'
        c, m = color_marks('slow-color')

        xs.append(x); ys.append(y)
        ls.append(lstyle); colors.append(c); markers.append(m); labels.append(l)

    color_marks.reset_counters(); color_marks.cind += 1
    for i in range(len(Mruns)):
        x = Mruns[i][:,ic]
        y = Mruns[i][:,it]

        lmax = int(Mres[i][0]); nr = int(Mres[i][1])
        l = r"MagIC, $\ell_\mathrm{{max}}={{{}}}$ $N_r={{{}}}$".format(lmax,nr)

        lstyle = '-'
        c, m = color_marks('slow-color')

        xs.append(x); ys.append(y)
        ls.append(lstyle); colors.append(c); markers.append(m); labels.append(l)

    xmin = np.min([np.min(i) for i in xs])
    xmax = np.max([np.max(i) for i in xs])
    xlim = (0.5*xmin, 200*xmax)

    MakePlot(xs, ys, labels, title, xlabel, ylabel, output, ls, colors, markers,
             dpi=dpi, legend=True, ylog=True, xlim=xlim, scaling_line=True)

if __name__ == "__main__":

    from docopt import docopt
    from glob import glob

    args = docopt(__doc__)

    Rfiles  = args['--r-files']
    Sfiles  = args['--s-files']
    Mfiles  = args['--m-files']
    output  = args['--output']
    dpi     = float(args['--dpi'])

    # expand any glob characters
    if (Rfiles is not None):
        Rfiles = glob(Rfiles)
    if (Sfiles is not None):
        Sfiles = glob(Sfiles)
    if (Mfiles is not None):
        Mfiles = glob(Mfiles)

    main(Rfiles, Sfiles, Mfiles, output, dpi)

