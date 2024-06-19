import numpy as np
import multiprocessing as mp
from events import Event

def powerlaw(logminval, logmaxval, power):

    """
    Draw a random number from a power law
    """

    logmin = logminval
    logmax = logmaxval

    c = -1.0 * logmax * power
    nmax = 10.0 ** (power * logmin + c)

    while True:
        log = np.random.uniform(logmin, logmax)
        n = 10.0 ** (power * log + c)

        if nmax * np.random.random() <= n:
            break

    return 10.0 ** log

def generate_events_worker(sub_events, freqs, peak_freq_pars, width_freq_pars, slope, peak_freq_dist, width_freq_dist):
    events = []
    for _ in range(sub_events):
        F = powerlaw(logFmin, logFmax, slope) * u.Jy * u.ms
        if peak_freq_dist == "uniform":
            f0 = np.random.uniform(peak_freq_pars[0].to(u.MHz).value, peak_freq_pars[1].to(u.MHz).value) * u.MHz
        elif peak_freq_dist == "normal":
            f0 = np.random.normal(peak_freq_pars[0].to(u.MHz).value, peak_freq_pars[1].to(u.MHz).value)  * u.MHz

        if width_freq_dist == "uniform":
            wf = np.random.uniform(width_freq_pars[0].to(u.MHz).value , width_freq_pars[1].to(u.MHz).value) * u.MHz
        elif width_freq_dist == "normal":
            wf = np.random.normal(width_freq_pars[0].to(u.MHz).value , width_freq_pars[1].to(u.MHz).value) * u.MHz

        event = Event(freqs=freqs, fluence=F, f0=f0, wf=wf)
        events.append(event)

    return events



def generate_events(freqs,
    peak_freq_pars,
    width_freq_pars,
    nevents = 100,
    logFmin = -2,
    logFmax = 3,
    slope = -1.5,
    peak_freq_dist = "uniform",
    width_freq_dist = "uniform",
    num_processes = 4):

    """
    Generate a list of toy FRB events
    """

    if num_processes is None:
        num_processes = mp.cpu_count()

    events = []
    pool = mp.Pool(processes=num_processes)

    # Calculate the number of events for each process
    events_per_process = nevents // num_processes
    remaining_events = nevents % num_processes

    # Generate events using multiprocessing
    results = [pool.apply_async(generate_events_worker, args=(events_per_process + (1 if i < remaining_events else 0), freqs, peak_freq_pars, width_freq_pars, slope, peak_freq_dist, width_freq_dist)) for i in range(num_processes)]
    pool.close()
    pool.join()

    # Collect results from processes
    for result in results:
        events.extend(result.get())

    return events


def make_subbands(nchan, a = 2, Z = 2, shifting = True ):

    """
    Make sub-bands like in (Trudu et al. 2024)
    """

    subbands = []

    for levelidx in range( Z + 1):

        nsubbands = int(a**(levelidx))
        chanpersub= int(a**(-levelidx)*nchan)

        #Regular sub-bands
        for chan in range(nsubbands):

            cstart = round((chan) * chanpersub)
            cstop  = round(cstart + chanpersub - 1)

            subbands.append([cstart, cstop])

        # Shifted sub-bands
        if shifting == True:
            for chan in range(nsubbands - 1):
                cstart = round((chan) * chanpersub  + chanpersub / 2)
                cstop  = round(cstart + chanpersub - 1)

                subbands.append([cstart, cstop])

    subbands = np.array(subbands)

    return subbands

def search(events,
          receiver,
          fstart = None,
          fend = None,
          threshold = 10,
          verbose = False):

    """
    Evaluate if a FRB will be detected in a given sub-band
    """

    bw    = receiver.bw
    fc    = receiver.fc
    Fmin  = receiver.get_mdf(threshold = threshold)
    freqs = receiver.get_freqs()

    if fstart is None:
        fstart = freqs[0]
    if fend is None:
        fend = freqs[-1]

    results = []


    for event in events:
        Fmeasured = event.fluence_perband(fstart, fend) * u.Jy * u.ms

        if Fmeasured.value >= Fmin.value * ( (fend.value - fstart.value) / bw.value ):
            results.append(True)
            if verbose == True:
                print("Burst Detected")
        else:
            results.append(False)
            if verbose == True:
                print("Burst Missed")
    return np.array(results)

def subbanded_search(events, receiver, subbands, verbose=False, threshold=10, batch_size=100):

    """
    Perform a sub-banded search in a given list of events
    """
    
    bw = receiver.bw
    fc = receiver.fc
    Fmin = receiver.get_mdf(threshold=threshold)
    freqs = receiver.get_freqs()
    nchan = receiver.nchan

    results = []
    event_batches = [events[i:i + batch_size] for i in range(0, len(events), batch_size)]

    def process_batch(batch):
        batch_results = []
        for event in batch:
            results_subband = []
            for row in subbands:
                cstart, cstop = row[0], row[1]
                fstart, fend = freqs[cstart], freqs[cstop]
                results_perband = search([event], receiver, fstart=fstart, fend=fend, threshold=threshold, verbose=False)
                results_subband.append(np.sum(results_perband))
            batch_results.append(any(results_subband))
        return batch_results

    pool = mp.Pool(mp.cpu_count())
    results = pool.map(process_batch, event_batches)
    pool.close()
    pool.join()

    results = [item for sublist in results for item in sublist]

    if verbose:
        for i, burst_detected in enumerate(results):
            if burst_detected:
                print(f"Burst Detected for event {i}")
            else:
                print(f"No Burst Detected for event {i}")

    return results
