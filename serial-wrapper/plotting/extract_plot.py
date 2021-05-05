"""
Extract performance data from a collection of files and plot it

Usage:
    extract.py [options]

Options:
    --r-files=<r>    Collection of pure Rayleigh files to process
    --s-files=<s>    Collection of Rayleigh-using-SHTns files to process
    --m-files=<m>    Collection of MagIC files to process
    --output=<o>     Set the output filename [default: scaling.png]
    --dpi=<d>        Resolution of image [default: 250]
    --single-nr=<n>  Only plot cases with the given Nr
    --method=<m>     Choose what to plot, best, mean, median [default: median]
"""

from __future__ import print_function
import matplotlib.pyplot as plt
plt.rcParams.update({"text.usetex":True})
import numpy as np
from plotting import MakePlot
from utils import ColorMarks
from glob import glob

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

def SplitIntoRuns(data, inds, method='best'):
    """
    Reorganize data into collection of (run vs core) entries
    data is (N_all_runs, Nentries), inds is dictionary of indices into data

    returns:
        -list of (Nrun, 3) arrays at a particular [lmax, nr], 0=Ncore, 1=time, 2=iter/sec
        -list of [lmax, Nr] items for each array
    """
    if (data is None): return [], []

    ir = inds['nr']; il = inds['lmax']; ic = inds['cores']; it = inds['time']; ii = inds['iter/sec']
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

    # choose best result for each core count & resolution
    for i, entry in enumerate(results):
        lmax = resolution[i][0]
        nr   = resolution[i][1]

        all_nc = list(set(entry[:,ic])) # unique/sorted list of core counts
        all_nc.sort()
        Nc = len(all_nc)

        # allocate storage
        result = np.zeros((Nc,3))

        for k, nc in enumerate(all_nc):
            K = np.where((entry[:,ic] == nc))[0] # pull out all entries at this lmax/nr/n_core

            vals = entry[K,:] # (Nrun, Nq)

            time  = vals[:,it] # all walltimes for this lmax/nr/n_core
            iters = vals[:,ii] # all iter/sec for this lmax/nr/n_core

            # choose best value for this lmax/nr/n_core
            if (method in ['best']):
                j = np.argmin(time)
                best_time = vals[j,it]

                j = np.argmax(iters)
                best_iter = vals[j,ii]

            elif (method in ['mean', 'avg']):
                best_time = np.mean(time)
                best_iter = np.mean(iters)

            elif (method in ['median']):
                best_time = np.median(time)
                best_iter = np.median(iters)

            else:
                raise ValueError("Unrecognized method = {}".format(method))

            result[k,0] = nc
            result[k,1] = best_time
            result[k,2] = best_iter

        # overwrite this entry
        results[i] = result

    return results, resolution

def main(Rfiles, Sfiles, Mfiles, output, dpi, single_nr, method):

    test = [x is None for x in [Rfiles,Sfiles,Mfiles]]
    if (all(test)):
        print("\nMust provide at least one of --r-files, --s-files, --m-files\n")
        return

    # set data and indices
    ind = {'lmax':0, 'nr':1, 'iter/sec':2, 'time':3, 'cores':4, 'nprow':5, 'npcol':6}
    print("\nParsing data files...")
    Rdata = ParseFiles(Rfiles) # (Nruns, Nentry)
    Sdata = ParseFiles(Sfiles) # (Nruns, Nentry)
    Mdata = ParseFiles(Mfiles, magic=True) # (Nruns, Nentry)

    # split into unique runs, i.e., same resolution but different core count
    print("Reorganizing data...")
    Rruns, Rres = SplitIntoRuns(Rdata, ind, method=method)
    Sruns, Sres = SplitIntoRuns(Sdata, ind, method=method)
    Mruns, Mres = SplitIntoRuns(Mdata, ind, method=method)

    # make plot
    xlabel = r"$N_\mathrm{cores}$"
    ylabel = "Walltime (sec)"
    title = "Scaling Results"

    xs = []; ys = []; labels = []; ls = []; colors = []; markers = []

    print("Adding data to plot...")
    color_marks = ColorMarks()
    ic = 0; it = 1
    for i in range(len(Rruns)):
        x = Rruns[i][:,ic]
        y = Rruns[i][:,it]

        lmax = int(Rres[i][0]); nr = int(Rres[i][1])
        l = r"Rayleigh, $\ell_\mathrm{{max}}={{{}}}$ $N_r={{{}}}$".format(lmax,nr)

        lstyle = '-'
        c, m = color_marks('slow-color')

        if ((single_nr is not None) and (not(np.allclose(nr, float(single_nr))))): continue

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

        if ((single_nr is not None) and (not(np.allclose(nr, float(single_nr))))): continue

        xs.append(x); ys.append(y)
        ls.append(lstyle); colors.append(c); markers.append(m); labels.append(l)

    color_marks.reset_counters(); color_marks.cind += 2
    for i in range(len(Mruns)):
        x = Mruns[i][:,ic]
        y = Mruns[i][:,it]

        lmax = int(Mres[i][0]); nr = int(Mres[i][1])
        l = r"MagIC, $\ell_\mathrm{{max}}={{{}}}$ $N_r={{{}}}$".format(lmax,nr)

        lstyle = '-'
        c, m = color_marks('slow-color')

        if ((single_nr is not None) and (not(np.allclose(nr, float(single_nr))))): continue

        xs.append(x); ys.append(y)
        ls.append(lstyle); colors.append(c); markers.append(m); labels.append(l)

    xmin = np.min([np.min(i) for i in xs])
    xmax = np.max([np.max(i) for i in xs])
    xlim = (0.5*xmin, 400*xmax)

    print("Drawing plot...\n")
    MakePlot(xs, ys, labels, title, xlabel, ylabel, output, ls, colors, markers,
             dpi=dpi, legend=True, ylog=True, xlim=xlim, scaling_line=True, width=6, height=None)

def Unglob(filenames):
    if ("," in filenames):
        files = []
        entries = filenames.split(",")
        for e in entries:
            files.extend(glob(e))
    else:
        files = glob(filenames)

    return files

if __name__ == "__main__":

    from docopt import docopt

    args = docopt(__doc__)

    Rfiles  = args['--r-files']
    Sfiles  = args['--s-files']
    Mfiles  = args['--m-files']
    output  = args['--output']
    dpi     = float(args['--dpi'])
    single_nr = args['--single-nr']
    method = args['--method']

    # expand any glob characters
    if (Rfiles is not None):
        Rfiles = Unglob(Rfiles)
    if (Sfiles is not None):
        Sfiles = Unglob(Sfiles)
    if (Mfiles is not None):
        Mfiles = Unglob(Mfiles)

    main(Rfiles, Sfiles, Mfiles, output, dpi, single_nr, method)

