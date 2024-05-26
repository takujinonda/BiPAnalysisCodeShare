import numpy as np
import pandas as pd

from scipy.signal import find_peaks, remez, filtfilt, spectrogram
from scipy.signal.windows import hamming
from scipy import stats

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
    Generate acceleration features (pitch, ODBA, and their 10 second moving means)

    Args:
        acc:                dataframe of acceleration signals.
        long_acc_name:      list of names for acceleration signals. Must conform to the following order - longitudinal (along body axis), dorsoventral (vertical axis), lateral.
        passb:              passband frequency (Hz).
        stopb:              stopband frequency (Hz). Must be higher than pass as lowpass filter is applied.
        fs:                 acceleartion signal sampling frequency (Hz).

    Returns:
        'out', Pandas Dataframe containing the following columns: 'pitch', pitch angles of acceleration signal, calculated via method from Sat et al. (2003). 'ODBA', Overall Dynamic Body Acceleration (the sum of absolute dynamic acceleration magnitudes). 'pitmn' and 'ODmn', 10 second moving means of pitch and ODBA respectively
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
    out['movsg'] = out.surge.rolling(2*fs, closed="both", min_periods=1).var()

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

def rollingSpecSum(spec,f,minFreq,maxFreq,fs=25,dur=60,inclusive=False):
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

def flightestimate(signal,fs,behav_data=None,dt=None,gap=30,cat='AT',low=3,high=5,removeErr=False,window=60,minGap=5,numPoints=20):
    """Estimate flight periods from dorsoventral acceleration. Flight is estimated by comparing frequencies within the `low` to `high` domain to `high`+. If `removeErr` is True, remove erroneous periods within a `gap`. The `numPoints` most likely `window` second periods of flight are extracted with at least `minGap` minutes inbetween.
    """
    f,_,Sxx = hammingSpect(signal,fs)
    rollSum = rollingSpecSum(Sxx, f, low, high, fs)
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
    if peaks[0] > troughs[0]:
        peaks = np.delete(peaks,0)
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

def flap(sig,fs,bout_gap=10,flap_freq=4,find_in_flight_periods=False,flinds=None,behav_data=None,dt=None,numPoints=10):
    """Find flapping signals in dorosventral signal `sig`. Flapping is extracted through peak-trough differences being greater than the inter-peak trough of the signal magnitude differences between maxima. These 'large' peaks and troughs are then grouped if they occur within half the typical flapping frequency `flap_freq`.

    Args:
        sig         - 
        fs          - 
        flap_freq   - 
        bout_gap    - 


    Returns:
    Two arrays of equal length as `sig`, one of flapping sequences `flap_mask` and another of flapping bouts `flap_bouts` (flap 1, not 0).
    """

    if find_in_flight_periods:
        _,_,flap_sig = peak_trough(sig)
        tst = flap_sig + flinds
        peaks = np.where(tst == 2)[0]
        troughs = np.where(tst == 4)[0]
    else:
        peaks,troughs,_ = peak_trough(sig)
    data = sig[peaks].values - sig[troughs].values
    ipt = interpeaktrough(data)
    
    # if only using estimated flight periods, recalculate peaks and troughs
    peaks,troughs,_ = peak_trough(sig)
    # retain only sufficiently large z displacements
    large_inds = (sig[peaks].values - sig[troughs].values) > ipt
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
    return flap_mask.astype(int), flap_bouts.astype(int)