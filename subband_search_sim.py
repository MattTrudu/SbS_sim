import os
import sys
import numpy as np
import astropy.units as u
import yaml
from receiver import Receiver
from utils import generate_events, subbanded_search


def main(configfile):

    with open(configfile, 'r') as yaml_file:
        config = yaml.safe_load(yaml_file)

    tsys  = config["tsys"] * u.K
    gain  = config["gain"] * u.K / u.Jy
    npol  = config["npol"]
    nchan = config["nchan"]
    fc    = config["fc"] * u.MHz
    bw    = config["bw"] * u.MHz
    dt    = config["dt"] * u.us
    name  = config["name"]

    nevents = config["nevents"]
    ncpu    = config["ncpu"]
    peak_freq_pars = config["peak_freq_pars"] * u.MHz
    width_freq_pars = config["width_freq_pars"] * u.MHz
    slope = config["slope"]
    peak_freq_dist = config["peak_freq_dist"]
    width_freq_dist = config["width_freq_dist"]
    logFmin = config["logFmin"]
    logFmax = config["logFmax"]

    a = config["a"]
    Z = config["Z"]
    shifting = config["shifting"]
    threshold = config["threshold"]
    batch_size = config["batch_size"]
    verbose = config["verbose"]
    outdir = config["outdir"]
    outname = config["outname"]

    if outdir == None:
        outdir = os.getcwd()

    freqs = np.linspace(peak_freq_pars[0] , peak_freqs_pars[1], 2 * nchan)

    receiver = Receiver(name = name, gain = gain, tsys = tsys, npol = npol, nchan = nchan, fc = fc, bw = bw, dt = dt)

    events = generate_events(freqs,peak_freq_pars,width_freq_pars,
                             nevents = nevents, logFmin = logFmin, logFmax = logFmax,
                             slope = slope, peak_freq_dist = peak_freq_dist,
                             width_freq_dist = width_freq_dist, num_processes = ncpu)

    subbands = make_subbands(nchan, a = a, Z = Z, shifting = shifting)

    results = subbanded_search(events, receiver, subbands, verbose = verbose, threshold = treshold, batch_size = batch_size)

    outname = os.path.join(outpath, outname)
    np.save(outname, results)

if __name__ == '__main__':

    configfile = sys.argv[1]

    main(configfile)
