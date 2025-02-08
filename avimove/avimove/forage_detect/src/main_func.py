#%%
import re
import numpy as np
import pandas as pd

from scipy import signal
from itertools import compress
from scipy.signal import find_peaks
from scipy import stats

def readAxy(filename, delim = "\t", cols = None, datetimeFormat = "%d/%m/%Y %H:%M:%S", ): 
    """
    Read in AxyTrek GPS data (txt files) as output by X Manager

    Args:

        filename:   path to AxyTrek txt file
        cols:       string to denote if acceleration and GPS ('acc') or GPS only ('gps') should be read in
        colnames:   list of string column names to be assigned. Must be of same length as cols. Defaults to ['Date', 'Time', 'lat', 'lon']

    Returns:
        Pandas dataframe of columns `colnames`. A formatted DateTime column (named DT) is generated.
    """

    accCols = ['Timestamp','X','Y','Z']
    gpsCols = ['Timestamp','location-lat','location-lon']

    # read in based on requested columns
    if cols.lower() == 'gps':
        df = pd.read_csv(filename, sep = ",", header = 0, usecols = gpsCols)
    elif cols.lower() == 'acc':
        df = pd.read_csv(filename, sep = ",", header = 0, usecols = accCols)
    else:
        df = pd.read_csv(filename, sep = ",", header = 0)
    df = pd.read_csv(filename, sep = delim, usecols = cols)
    df.DT = [dtFormat(x) for x in df.DT] # ensure correct datetime formats
    df['DT'] = pd.to_datetime(df['Date'] + " " + df['Time'], format = datetimeFormat)
    return df

def readDVL(filename,accStart=None,vidStart=None,removePriorVid=False):
    """
    Read Little Leonardo DVL data logger (400M) data. File locations and start-times of acceleration and video recording (dd/mm/yyyy HH:MM:SS) required.
    """
    # read in acceleration data
    acc = pd.read_table(filename, skiprows=7, sep=',', usecols=[0,1,2])
    # remove header whitespace
    acc.rename(columns=lambda x: x.strip(), inplace=True)
    acc['DT'] = pd.date_range(accStart, periods=len(acc), freq=f"{1/self.fs * 10**3}ms")# add time series
    if removePriorVid:
        # select data within video range
        acc = acc[(acc.DT >= vidStart) & (acc.DT < (vidStart + pd.Timedelta(hours=2)))]

def readBIP(filename,cols=None):
    """
    Read in BiP-formatted data

    Args:

        filename:   string of full path to file
        cols:       string value indicating if number of columns to read in can be reduced. For GPS and acceleration, use 'acc', for GPS only, use 'gps'. Defaults to None.

    Returns:
        Pandas dataframe of BIP system formatted data. Columns returned depend on `cols` argument. 'acc' returns formatted datetime (DT), latitude (lat), longitude (lon), longitudinal acceleration (X), lateral acceleration (Y), and dorsoventral acceleration (Z). 'gps' returns datetime, latitude, and longitude. If no cols argument given, all data is read in, including pressure (pressure), temperature (temp), height about sea level (altitude), and ground speed (spd).
    """

    accCols = ['time','latitude','longitude','acceleration_longitudinal','acceleration_lateral','acceleration_dorso_ventral']
    gpsCols = ['time','latitude','longitude']

    # read in based on requested columns
    if cols.lower() == 'gps':
        df = pd.read_csv(filename, sep = ",", header = 0, usecols = gpsCols)
    elif cols.lower() == 'acc':
        df = pd.read_csv(filename, sep = ",", header = 0, usecols = accCols)
    else:
        df = pd.read_csv(filename, sep = ",", header = 0)
         
    # df = pd.read_csv(filename, sep = ",", header = 0, usecols = cols)
    # rename columns for later use
    df.rename(columns = {'time':'DT','latitude':'lat','longitude':'lon','acceleration_longitudinal':'X','acceleration_lateral':'Y','acceleration_dorso_ventral':'Z','pressure':'pressure','temperature':'temp','height_above_mean_sea_level':'altitude','ground_speed':'spd'}, inplace = True)
    
    df.DT = [dtFormat(x) for x in df.DT] # ensure correct datetime formats
    df['DT'] = pd.to_datetime(df['DT'], format = "%Y-%m-%d %H:%M:%S.%f")

    return df

def gps_distanceSingle(longitudes, latitudes, latVal, lonVal):                                                                       
    # taken from https://www.tjansson.dk/2021/03/vectorized-gps-distance-speed-calculation-for-pandas/, thanks to Thomas Jansson for use of this function
    """                                                                                                                 
    Calculates distances between an array of GPS points and a single position via a vectorized haversine calculation the great circle distance between two arrays of points on the earth (specified in decimal degrees). All args must be of equal length.                                         
 
    Args:                                                                                                               
        longitudes: pandas series of longitudes                                                                         
        latitudes:  pandas series of latitudes                               
        latVal:     latitude value of single position                                           
        lonVal:     longitude value of single position
 
    Returns:                                                                                                            
        Distance is returned an array in meters.                                                                             
 
    Example:                                                                                                            
        >>> df['gpsDist'] = gps_distance(df.longitude, df.latitude, 39.400,141.998)
    """
 
    # lon1 = longitudes[:-1].reset_index(drop = True)                                                                                       
    # lat1 = latitudes[:-1].reset_index(drop = True)                                                                                 
    # lon2 = longitudes[1:].reset_index(drop = True)                                                                                       
    # lat2 = latitudes[1:].reset_index(drop = True)                                                                                 

    return [distance.distance((x,y),(latVal,lonVal)).km for x,y in zip(latitudes,longitudes)]

# to be placed elsewhere
import re

def dtFormat(x):
    """
    Format incoming datetimes from BiP system. Typically datetimes are brought in the format YYYY-mm-dd HH:MM:SS+00:00, however, SS can contain decimal values or not. This function returns the same datetime in string format with a decimal place added if none is originally present and removes '+00:' from the string so that datetimes can be converted to datetime format in Pandas.

    Args:
    
        x:  datetime in string format

    Returns:
        String of x in correct datetime format %Y-%m-%d %H:%M:%S.%f
    """

    # extract all before '+'
    out = re.findall('^(.*)\+',x)[0]
    # remove all non-datetime characters
    out = re.sub('[^0-9 .:-]','',out)
    if '.' not in out: # add milliseconds if not present (3 d.p.)
        out = out + '.000'
    if out.count('.') > 1:
        out = re.findall('[^.]*.[^.]*',out)[0] # extract all before second period
    # ensure 3 d.p.
    msLength = len(re.findall('[^.]*.',out)[1])
    if msLength != 3:
        if msLength < 3:
            out = out + ('0' * (3-msLength))
        else:
            out = out[:-(msLength - 3)]

    return out

def gps_speed(longitudes, latitudes, timestamps):                                                                       
    # taken from https://www.tjansson.dk/2021/03/vectorized-gps-distance-speed-calculation-for-pandas/, thanks to Thomas Jansson for use of this function
    """                                                                                                                 
    Calculates the instantaneous speed from the GPS positions and timestamps. The distances between the points          
    are calculated using a vectorized haversine calculation the great circle distance between two arrays of points on   
    the earth (specified in decimal degrees). All args must be of equal length.                                         
 
    Args:                                                                                                               
        longitudes: pandas series of longitudes.                                                                         
        latitudes:  pandas series of latitudes.                                                                          
        timestamps: pandas series of timestamps.                                                                         
 
    Returns:                                                                                                            
        Speed is returned an array in m/s.                                                                             
 
    Example:                                                                                                            
        >>> df['gpsSpeed'] = gps_speed(df.longitude, df.latitude, df.recordedAt).
    """
 
    lon1 = longitudes[:-1].reset_index(drop = True)                                                                                       
    lat1 = latitudes[:-1].reset_index(drop = True)                                                                                 
    lon2 = longitudes[1:].reset_index(drop = True)                                                                                       
    lat2 = latitudes[1:].reset_index(drop = True)                                                                                 
 
    # Vectorized haversine calculation                                                                                  
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])                                                  
    a = np.sin((lat2 - lat1) / 2.0)**2 + (np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2.0)**2)                 
    dist = np.array(6371 * 2 * np.arcsin(np.sqrt(a)) * 1000)
    time_array = (timestamps.diff().dt.seconds ).values[1:]                                                       
 
    # Calculate the speed                                                                                               
    time_array[time_array == 0] = np.nan  # To avoid division by zero                                                   
    speed = np.array(dist / time_array)

    # Make the arrays as long as the input arrays                                                                        
    speed = np.insert(speed, 0, np.nan, axis=0)
    dist = np.insert(dist, 0, np.nan, axis=0)                                                               
    return dist,speed

def distSpeed(lat,lon,DT,threshold=None):
    """
    Calculates distances (metres) and speed (m/s) of GPS positions and datetime information. Uses a speed threshold (default None) to define erroneous GPS positions. GPS and datetime lengths must agree.

    Args:

        lat:    array of latitude points in decimal degrees.
        lon:    array of longitude points in decimal degrees.
        DT:     datetime array in Pandas datetime format.

    Returns:
        Float arrays of distance and speed values.
    """

    dist,speed = gps_speed(lat,lon,DT)
    if threshold is not None:

        while np.nanmax(speed) > threshold:

            # remove erroneous GPS values
            lon[speed > threshold] = np.nan
            lat[speed > threshold] = np.nan

            # recalculate speed
            dist,speed = gps_speed(lat,lon,DT)
    
    return dist,speed

# def main_func(df):
#     df['time'] = pd.to_datetime(df['time'], format='ISO8601')
#     df = df.drop_duplicates('time', keep='last')
#     df = df.dropna()
#     df = df.reset_index()
#     df = df.drop(['index'], axis=1)
#     df.rename(columns={'time': 'DT', 'latitude': 'lat', 'longitude': 'lon'}, inplace=True)
#     df['DT'] = [dtFormat(str(x)) for x in df['DT']]
#     df['DT'] = pd.to_datetime(df['DT'],format="%Y-%m-%d %H:%M:%S.%f") #一旦、+00:00を除くために文字列にする
    
#     if len(df) >0:
#         out = windEstimation2(df, isBp = True)
#         out = out.rename(columns={"Time":"time", "Lat":"latitude", "Lon":"longitude"})
        
#         #すでにUTCなので不要
#         out['time'] = pd.to_datetime(out['time'], format='%Y-%m-%d %H:%M:%S.%f', errors='coerce').dt.tz_localize('UTC')

#         return out
#     else:
#         return pd.DataFrame(columns=['Time','Lat','Lon','MeanHead','X','Y'])

def removeNearCS(data, captureSite, distThreshold = 1.5):
    """
    Remove acceleration and GPS data captured within a set distance of the capture site.

    Args:
        data:           pandas dataframe of tag data as formatted by readBIP. Expected columns are DT (datetime, formatted), lat (latitude in decimal degrees), and lon (longitude in decimal degrees).
        captureSite:    capture site latitude and longitude in decimal degrees.
        distThreshold:  desired distance threshold from capture site in km. Defaults to 1.5.

    Returns:
        Pandas dataframe of similar format as data with all points located within distThreshold of captureSite removed.
    """
    # highest data sample rate (Hz)
    fs = int(np.timedelta64(1,'s') / stats.mode(np.diff(data.DT),keepdims=False)[0])

    # GPS sample rate (fixes per minute)
    gpsFS = 60/int(np.timedelta64(1,'m') / stats.mode(np.diff(data.dropna()['DT']),keepdims=False)[0])

    # find distances from capture site
    capSiteDist = np.array(gps_distanceSingle(data.lon.dropna(),data.lat.dropna(),captureSite[0],captureSite[1]))

    inds = data.lat.dropna().index[np.where(capSiteDist < distThreshold)[0]] # find indeces to be removed from original dataset

    # if any difference in GPS and acc sampling rates means we must extend GPS indeces
    samplerDiff = int(np.floor(fs * gpsFS / 2)) # number of fs samples either side of GPS position

    # create ranges to be removed
    rem = []
    for x in inds:
        rem.extend(np.arange(x - samplerDiff, x + samplerDiff))
    rem = np.array(rem)

    return data.drop(rem)

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
    eqFil=signal.remez(101,[0,passband,stopband,fs*.5],[1,0],fs=fs)
    # return 'static' signals for each provided signal
    static = signal.filtfilt(b=eqFil,a=1,x=sig)
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
    win = signal.windows.hamming(winsize);

    f, t, Sxx = signal.spectrogram(sig,fs,window=win,noverlap=numoverlap)

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

def flightestimate(signal,fs,behav_data=None,dt=None,gap=30,cat='AT',low=3,high=5,removeErr=True,window=60,minGap=5,numPoints=20):
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

def flap(sig,fs,bout_gap=10,flap_freq=4,find_in_flight_periods=False,behav_data=None,dt=None,numPoints=10):
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
        _, flinds = flightestimate(sig,20,behav_data=behav_data,dt=dt,numPoints=numPoints)
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

def flight_pitch_changes(sig,fl_inds):
    _, _, pit_sig = peak_trough(sig)
    # find overlap with flights
    fl_pit_pk = np.where(pit_sig + fl_inds == 2)
    fl_pit_trgh = np.where(pit_sig + fl_inds == 2)
# %%