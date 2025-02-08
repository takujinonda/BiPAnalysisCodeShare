import numpy as np
import pandas as pd

from scipy.signal import find_peaks, remez, filtfilt, spectrogram
from scipy.signal.windows import hamming
from scipy import stats
from itertools import compress
from statistics import median

def lowEquiFilt(sig, passband, stopband, fs):
    """
    Generate and apply lowpass equiripple filter (equivalent to MATLAB default lowpass filter) to produce static (low pass) and dynamic (original - low pass - see Patterson et al. 2019)

    Args:
        sig:        signal to apply filter to.
        passband:   passband (Hz).
        stopband:   stopband (Hz).
        fs:         sig sampling frequency (Hz).

    Returns:
        `static`, the low pass filtered signal, and `dynamic`, the difference between the original signal and `static`. Both are arrays of same length as sig
    """
    # generate equiripple filter
    eqFil= remez(101,[0,passband,stopband,fs*.5],[1,0],fs=fs)
    # return 'static' signals for each provided signal
    static = filtfilt(b=eqFil,a=1,x=sig)
    # return 'dynamic' acceleration signal
    dynamic = np.array(sig - static)
    
    return static, dynamic

# create static and dynamic (1 pass, 1.5 stop), pitch, ODBA
def accFeatures(acc, long_acc_name, passb, stopb, fs):
    """
    Generate acceleration features (pitch, ODBA, and their 10 second moving
    means)

    Args
    ----
    acc
        dataframe of acceleration signals.
    long_acc_name
        list of names for acceleration signals. Must conform to the following
        order - longitudinal (along body axis), dorsoventral (vertical axis),
        lateral.
    passb
        passband frequency (Hz).
    stopb
        stopband frequency (Hz). Must be higher than pass as lowpass filter is
        applied.
    fs
        acceleartion signal sampling frequency (Hz).

    Returns:
        'out', Pandas Dataframe containing the following columns: 'pitch', pitch
        angles of acceleration signal, calculated via method from Sat et al.
        (2003). 'ODBA', Overall Dynamic Body Acceleration (the sum of absolute
        dynamic acceleration magnitudes). 'pitmn' and 'ODmn', 10 second moving
        means of pitch and ODBA respectively
    """

    # take acceleration signal names
    sLong, dLong = lowEquiFilt(acc[long_acc_name[0]], passb, stopb, fs)
    _, dZ = lowEquiFilt(acc[long_acc_name[1]], passb, stopb, fs)
    _, dLat = lowEquiFilt(acc[long_acc_name[2]], passb, stopb, fs)

    pitch = np.arcsin(np.clip(sLong,-1,1)) * 180/np.pi
    Odba = sum([np.abs(dLong),np.abs(dZ),np.abs(dLat)])

    out = pd.DataFrame({'pitch':pitch,'ODBA':Odba,'surge':dLong,'DZ':dZ})
    out['pitmn'] = out.pitch.rolling(10*fs, closed="both", min_periods=1).mean()
    out['ODmn'] = out.ODBA.rolling(10*fs, closed="both", min_periods=1).mean()
    # out['movsg'] = out.surge.rolling(2*fs, closed="both", min_periods=1).var()

    return out

def hammingSpect(sig,fs=25):
    """
    Generate spectrogram data of sig using a Hamming window of 4 seconds with 85% overlap.
    
    Args:
        sig:    signal to generate spectrogram.
        fs:     sig sampling frequency (Hz). Defaults to 25.
        
    """
    # generate spectrogram (4 second window, overlap of 85%, hamming window)

    #set window size of 4 seconds
    winsize = fs*4
    #set overlap between windows (will determine the temporal resolution
    numoverlap = np.floor(.9875*winsize).astype(int); #(85%)
    win = hamming(winsize);

    f, t, Sxx = spectrogram(sig,fs,window=win,noverlap=numoverlap)

    return f, t, Sxx

def rollingSpecSum(sig,minFreq,maxFreq,fs=25,dur=60,inclusive=False):
    """
    Create a rolling sum of the difference between spectrogram intensities from minFreq to maxFreq and intensities above maxFreq, , i.e. summed values between 3 and 5 Hz - summed values above 5 Hz.
    
    Args:
        spec:       spectrogram object.
        f:          array of sample frequencies.
        minFreq:    minimum frequency for sum (Hz).
        maxFreq:    maximum frequency for sum (Hz).
        fs:         spectrogram sampling frequency (Hz). Defaults to 25.
        dur:        rolling sum duration (s).
        inclusive:  should boundaries be included. Defaults to False.

    Returns:
        Summed frequency intensities across spec object.
    """
    f,_,spec = hammingSpect(sig,fs)
    if inclusive:
        spectDiff = pd.Series(np.sum(spec[(f >= minFreq) & (f <= maxFreq),:],axis=0) - np.sum(spec[f > maxFreq,:],axis=0))
    else:
        spectDiff = pd.Series(np.sum(spec[(f > minFreq) & (f < maxFreq),:],axis=0) - np.sum(spec[f > maxFreq,:],axis=0))
    return spectDiff.rolling(dur*fs, closed = "both", min_periods = 1).sum()

def maxWithGap(sig,fs,window=60,minGap=5,numPoints = 20):
    """
    Calculate the highest values across time-series signal `sig` with a minimum time gap between
    
    Args:
        sig:        signal whose largest values are to be found.
        fs:         signal sampling frequency (Hz).
        windows:    duration of moving window over which to search (mins). Defaults to 1.
        minGap:     minimum time gap between largest values (mins). Defaults to 5.
        numPoints:  number of 'max points' to be found. Defaults to 20.

    Returns:
        List of length numPoints of ranges indicating indeces of highest values within sig.        
    """
    out = []
    sigad = sig.copy()
    while len(out) != numPoints:
        # create index range around highest value
        out.append(np.arange(np.argmax(sigad) - round(fs*window/2) + fs*2,
                  np.argmax(sigad) + round(fs*window/2)) + fs*2)
        # reduce magnitude of this period and 5 minutes surrounding
        sigad[np.arange(np.max([0,np.argmax(sigad) - round(fs*window/2) - (fs*60*minGap)]),np.min([len(sigad),np.argmax(sigad) + round(fs*window/2) + (fs*60*minGap)]))] = np.min(sigad)
    return out

def reduceErroneous(signal,behav_data,dt,fs,gap=30,cat='AT'):
    """Reduce the magnitude of erroneous periods of data `signal` according to behaviour data `behav_data`. Surrounding `gap` seconds are reduce of the time period of concern.

    Args:
        signal      - Signal for which erroneous periods should have their magnitude reduced
        behav       - Behavioural data with classification `cat` for erroneous data
        dt          - datetime array of signal
        fs          - `signal` sampling frequency
        gap         - window over which to reduce magnitude (in seconds)
        cat         - behav category to indicate erroneous periods, defaults to 'AT'

    Returns:
    `out`, array of same size as `signal`, now with segments of reduced magnitude where `behav_data` indicated erroneous behaviour
    """
    sigad = signal.copy()
    if not behav_data.Behaviour.eq(cat).any():
        print("No erroneous data found")
    else:
        behAT = np.any([(dt >= (x - pd.Timedelta(gap,'sec'))) & (dt <= (x + pd.Timedelta(gap,'sec'))) for x in behav_data.Time[behav_data.index[behav_data.Behaviour == cat]].round("s").values], axis = 0)
        sigad.loc[behAT[(2*fs):-(2*fs)+1]] = min(sigad)
    return sigad

def flightestimate(signal,rollSum,fs,behav_data=None,dt=None,gap=30,cat='AT',removeErr=False,window=60,minGap=5,numPoints=20):
    """Estimate flight periods from dorsoventral acceleration. Flight is estimated by comparing frequencies within the `low` to `high` domain to `high`+. If `removeErr` is True, remove erroneous periods within a `gap`. The `numPoints` most likely `window` second periods of flight are extracted with at least `minGap` minutes inbetween.
    """
    if removeErr:
        rollSum = reduceErroneous(rollSum,behav_data,dt,fs,gap,cat)
    out = maxWithGap(rollSum,fs,window,minGap,numPoints)
    flight = np.zeros(len(signal))
    for x in out:
        flight[x] = 1
    return out, flight.astype(int)

def find_gaps(signal, gap_size):
    """Identify gaps in bool array `signal` greater than `gap_size` in length. Extend True segments to contain all elements contained within gap.

    Returns arrays of start and end points.
    """
    past_threshold = np.array([x-y for x, y in zip(signal[1:],signal)]) > gap_size
    past_threshold = np.array(list(compress(range(len(past_threshold)),past_threshold)))
    # ends = past_threshold.append(len(flapInds))
    ends = np.append(past_threshold,(len(signal)-1))
    starts = np.insert(past_threshold + 1,0,0)
    # mask = np.zeros(ends[-1]).astype(int)
    # for x,y in zip(starts,ends):
    #     mask[x:y] = 1
    return signal[starts],signal[ends]

def peak_trough(sig):
    """Extract peaks and troughs of signal `sig`. Will start with peak and end in trough to ensure equal size
    """
    peaks,_ = find_peaks(sig)
    troughs,_ = find_peaks(-sig)
    while len(peaks) != len(troughs):
        if peaks[0] > troughs[0]:
            troughs = np.delete(troughs,0)
        if peaks[-1] > troughs[-1]:
            peaks = np.delete(peaks,-1)

    # create alternative of bool array indicating presence/absence of peaks (1) and troughs (2)
    flap_sig = np.zeros(len(sig))
    flap_sig[peaks] = 1
    flap_sig[troughs] = 3
    return peaks, troughs, flap_sig

def interpeaktrough(mags):
    """Find the inter-peak trough of signal `mags`. Only the first two peaks are considered
    """
    kde = stats.gaussian_kde(mags)
    x = np.linspace(mags.min(),mags.max(),100)
    p = kde(x)

    ppeaks,_ = find_peaks(p)
    pks,_ = find_peaks(-p[ppeaks[0]:ppeaks[1]])

    return x[pks]

def flatten(xss):
    return [x for xs in xss for x in xs]

def peak_trough_in_flight(sig,fl_inds):
    """Identify peaks/troughs of signal `sig` within defined flight periods `fl_inds`
    """

    peaks,troughs,_ = peak_trough(sig)
    flpeaks = []
    fltroughs = []
    for inds in fl_inds:
        pks = np.sort(list(set(peaks).intersection(inds)))
        trghs = np.sort(list(set(troughs).intersection(inds)))
        if pks[0] > trghs[0]:
            trghs = np.delete(trghs,0)
        if pks[-1] > trghs[-1]:
            pks = np.delete(pks,-1)
        flpeaks.append(pks)
        fltroughs.append(trghs)
    peaks = flatten(flpeaks)
    troughs = flatten(fltroughs)
    return peaks, troughs

def flap(sig,fs,bout_gap=30,flap_freq=4,find_in_flight_periods=False,flinds=None):
    """
    Find flapping signals in dorosventral signal `sig`. Flapping is extracted through peak-trough differences being greater than the inter-peak trough of the signal magnitude differences between maxima. These 'large' peaks and troughs are then grouped if they occur within half the typical flapping frequency `flap_freq`.

    Parameters
    ----------
    sig
        Dorsoventral acceleration signal.
    fs
        Sampling frequency (Hz).
    flap_freq
        Typical expected frequency of flapping (Hz).
    bout_gap
        Maximum gap between segments to separate from 'glides' and group flapping periods into bouts.


    Returns:
    --------
    flap_mask
        Periods of flapping grouped by half the typical flapping frequency.
    flap_bouts
        Periods of flapping grouped by the bout_gap period.
    starts
        Start indeces of flapping bouts.
    ends
        End indeces of flapping bouts.
    """

    if find_in_flight_periods:
        peaks,troughs = peak_trough_in_flight(sig,flinds)
    else:
        peaks,troughs,_ = peak_trough(sig)
    data = np.absolute(sig[peaks].values - sig[troughs].values)
    ipt = interpeaktrough(data)
    
    # if only using estimated flight periods, recalculate peaks and troughs
    peaks,troughs,_ = peak_trough(sig)
    # retain only sufficiently large z displacements
    large_inds = np.absolute(sig[peaks].values - sig[troughs].values) > ipt
    # convert to flapping indeces (of acc signal)
    flap_inds = np.sort(np.concatenate([peaks[large_inds],troughs[large_inds]]))
    # find gaps between flap signals greater than twice the typical flapping frequency
    flap_mask = np.zeros(len(sig))
    start,end = find_gaps(flap_inds,2*fs/flap_freq)
    for x,y in zip(start,end):
        flap_mask[x:y] = 1
    
    flap_inds = np.array([i for i, x in enumerate(flap_mask) if x == 1])
    starts,ends = find_gaps(flap_inds, fs*bout_gap)
    flap_bouts = np.zeros(len(sig))
    for x,y in zip(starts,ends):
        flap_bouts[x:y] = 1
    return flap_mask.astype(int), flap_bouts.astype(int), starts, ends

def flight_est_thresholds(acc,fl_inds):
    """Define threshold values required for behaviour discrimination. Usage requires full acceleration data `acc` alongside indeces of estiamted flight periods `fl_inds`
    """
    pitVar = median([np.var(acc.pitch[x]) for x in fl_inds])
    pitFLmn = median([np.max(acc.pitMn[x]) for x in fl_inds])
    ODmFL = median([np.min(acc.ODmn[x]) for x in fl_inds])
    return pd.DataFrame({'pitVar': pitVar,'pitFlmn':pitFLmn, 'ODmFL': ODmFL})

def flight_pitch_changes(sig,fl_inds,findVal=None):
    outs = []
    for x in fl_inds:
        pk,trgh,_ = peak_trough(sig[x])
        pk = x[pk]
        trgh = x[trgh]
        if findVal == 'max':
            outs.append(np.max(np.array(sig[pk]) - np.array(sig[trgh])))
        elif findVal == 'min':
            outs.append(np.min(np.array(sig[pk]) - np.array(sig[trgh])))
    return median(outs)