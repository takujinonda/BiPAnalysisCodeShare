import cmath
import math
import re
from functools import partial
from statistics import mode

#import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sp
from scipy import stats
from scipy.optimize import minimize


def Likelihoodww(data1: np.ndarray,data2: np.ndarray,cv: np.ndarray):
    """                                                                                                                 
    Calculates log-likelihood function for testing model output                                         
 
    Args:                                                                                                               
        data1:  numpy array of speed values
        data2:  numpy array of headings
        cv:     estimated constant airspeed
 
    Returns:                                                                                                            
        log-likelihood function
 
    Example:                                                                                                            
        >>> minimize(Likelihoodww(speed,headings,constantV),initPars)
    """           

    def f(par):
        b = cv/sp.special.gamma(1+1/par[0])
        L = 0
        for i,j in zip(data1,data2):
            r1 = np.sqrt(pow((i*np.cos(j) - par[3]),2) + pow((i*np.sin(j) - par[4]),2))
            rx = (i*np.cos(j)-par[3])/r1
            ry = (i*np.sin(j)-par[4])/r1
            lp = (par[0]-2) * math.log(r1) - (r1/b)**par[0] + par[1]*rx + par[2]*ry + math.log(par[0]) - math.log(b) + (1-par[0])*math.log(b) - math.log(sp.special.iv(0,np.sqrt(par[1]**2 + par[2]**2),))
            L = L+lp
        return L/(-1)
    return f

def Weibull_sd(a,b):
    """Calculate the standard deviation of a Weibull distribution
    
    Args:

        a:  mean heading vector length
        b:  the shape parameter

    Returns:
        numpy float

    Example:
        >>> Weibull_sd(14,1.3)
    """

    return b*np.sqrt(sp.special.gamma(1+2/a) - sp.special.gamma(1+1/a)*sp.special.gamma(1+1/a))

def Weibull_mean(a,b):
    """Mean of Weibull distribution

    Args:

        a:  mean heading vector direction
        b:  shape parameter

    Returns:
        numpy float
    
    
    Example:
        >>> Weibull_mean(14,1.3)
    """
    return b*sp.special.gamma(1+1/a)

def Von_Mises_sd(kappa):
    """Standard deviation of von Mises distribution

    Args:
        kappa:  concentration parameter

    Returns:
        NumPy float

    Example:
        >>> Von_Mises_sd(1.4)
    """

    return 1/np.sqrt(kappa)

def readAxyGPS(filename, delim = "\t", cols = [0,1,2,3], colnames = ['Date', 'Time', 'lat', 'lon'], datetimeFormat = "%d/%m/%Y %H:%M:%S"): 
    """
    Read in AxyTrek GPS data (txt files) as output by X Manager

    Args:

        filename:   path to AxyTrek txt file
        cols:       list of columns to be read in, defaults to [0,1,2,3]
        colnames:   list of string column names to be assigned. Must be of same length as cols. Defaults to ['Date',    'Time', 'lat', 'lon']

    Returns:
        Pandas dataframe of columns `colnames`. A formatted DateTime column (named DT) is generated.
    """
    df = pd.read_csv(filename, sep = delim, usecols = cols,
    names = colnames)
    df.DT = [dtFormat(x) for x in df.DT] # ensure correct datetime formats
    df['DT'] = pd.to_datetime(df['Date'] + " " + df['Time'], format = datetimeFormat)
    return df

def readBIPAxy(filename):
    """
    Read in AxyTrek data as formatted by the BIP system

    Args:

        filename:   full path to file

    Returns:
        Pandas dataframe of all columns from the BIP system (datetime 'DT', latitude 'lat', and longitude 'lon'). A formatted DateTime column (named DT) is generated)
    """
    
    df = pd.read_csv(filename, sep = ",", header = 0).dropna().reset_index()
    df = df[['time','latitude','longitude']].rename(columns={'time':'DT','latitude':'lat','longitude':'lon'})
    df.DT = [dtFormat(x) for x in df.DT] # ensure correct datetime formats
    df['DT'] = pd.to_datetime(df['DT'], format = "%Y-%m-%d %H:%M:%S.%f")
    return df

def nearest(items, pivot):
    """Find the nearest value within an array

    Args:

        items:  array to be searched within
        pivot:  value to find the nearest within array items

    Returns:
        Member of items closest in value to pivot

    Example:
        >>> data = np.arange(1,100,3.5) # generate regular sequence with step of 3.5
        >>> nearValue = nearest(data, 27) # find the nearest value to 27
    """
    return min(items, key=lambda x: abs(x - pivot))

def nearestInd(items, pivot):
    """Find the nearest index within an array

    Args:
    
        items:          array to be searched within
        pivot:          value to find the index of the nearest value within items

    Returns:
        Index value for nearest member of items in value to pivot
    
    Example:
        >>> data = np.arange(1,100,3.5) # generate regular sequence with step of 3.5
        >>> nearIndex = nearestInd(data, 27) # find index of value within data closest to 27
    """

    return min(range(len(items)), key=lambda i: abs(items[i] - pivot))

def timeRescale(dat,tdiff,units='min'):
    """    
    Subset dataframe to reflect desired regular sampling interval. Nearest time values are used, `dat` must have datetime column 'DT'
    Args:

        dat:            pandas dataframe with datetime column'DT'  correctly formatted
        tdiff:          desired regular sampling interval
        units:          units of desired sampling intervals using pandas time/date conventions. Defaults to 'T'.

    Returns:
        Pandas dataframe resampled to desired regular sampling interval
    """

    # change indexing to datetime and resample
    out = dat.set_index('DT', drop = False).resample(str(tdiff) + units).nearest().dropna()
    # return resampled data with original times and standard indexing
    return out.set_index(pd.Index(range(len(out))))

def angles(longitudes,latitudes):
    """
    Angular differences between subsequent in latitude and longitude arrays

    Args:

        longitudes:     narray of longitudes in decimal degrees
        latitudes:      array of latitudes in decimal degrees

    Returns:
        Array of length len(longitudes) - 1 showing angles between consecutive GPS positions in radians

    Example:
        >>> long = np.arange(110,120,0.5)
        >>> lat = np.arange(38,39,0.05)
        >>> angularChange = angles(long,lat)
    """

    lon1 = longitudes[:-1].reset_index(drop = True)                                                                                       
    lat1 = latitudes[:-1].reset_index(drop = True)                                                                                        
    lon2 = longitudes[1:].reset_index(drop = True)                                                                                        
    lat2 = latitudes[1:].reset_index(drop = True)
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])   

    X = lon2 - lon1 
    Y = lat2 - lat1
    return np.append(np.nan,np.arctan2(Y,X))

def gps_speed(longitudes, latitudes, timestamps):                                                                       
    # taken from https://www.tjansson.dk/2021/03/vectorized-gps-distance-speed-calculation-for-pandas/, thanks to Thomas Jansson for use of this function
    """                                                                                                                 
    Calculates the instantaneous speed from the GPS positions and timestamps. The distances between the points          
    are calculated using a vectorized haversine calculation the great circle distance between two arrays of points on   
    the earth (specified in decimal degrees). All args must be of equal length.                                         
 
    Args:                                                                                                               
        longitudes: pandas series of longitudes                                                                         
        latitudes:  pandas series of latitudes                                                                          
        timestamps: pandas series of timestamps                                                                         
 
    Returns:                                                                                                            
        Speed is returned an array in m/s.                                                                             
 
    Example:                                                                                                            
        >>> df['gpsSpeed'] = gps_speed(df.longitude, df.latitude, df.recordedAt)
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

def dtFormat(x):
    """
    Format incoming datetimes from BiP system. Typically datetimes are brought in the format YYYY-mm-dd HH:MM:SS+00:00, however, SS can contain decimal values or not. This function returns the same datetime in string format with a decimal place added if none is originally present and removes '+00:' from the string so that datetimes can be converted to datetime format in Pandas.

    Args:
    
        x:  datetime in string format

    Returns:
        String of x in correct datetime foramt %Y-%m-%d %H:%M:%S.%f
    """
    # first, test if decimal place already present
    if bool(re.search('[.]',x)):
        x = re.sub('[+]','',x)
    else:
        x = re.sub('[+]','.',x)
    x = re.sub(':00$','',x)

    return x

def prePare(filename, convertToMin: bool = True, tdiff = 1, units = 'min', isBip: bool = True):
    """Prepare AxyTrek GPS data as per required for Goto original method. DateTime, distance, track speed and direction is added

    Args:
        filename:       full path for file to be read in
        convertToMin:   boolean for whether data should be resampled to regular time intervals
        tdiff:          resample time interval. Defaults to 1
        units:          resample time interval unit. Defaults to 'min'
        isBip:          boolean stating whether data is in BiP format

    Returns:
        Pandas dataframe of read in data. Adds columns 'dt' (elapsed time from fix to previous time point in seconds), 'dist' (distance travelled from previous point in m), 'track_speed' (in m/sec), 'track_direction' (direction between consecutive GPS positions in rad)
    """

    if isBip:
        df = readBIPAxy(filename)
    else:
        df = readAxyGPS(filename)
    if convertToMin:
        df = timeRescale(df, tdiff, units)
    df['dt'] = np.append(np.nan,(np.diff(df['DT']) / np.timedelta64(1,'s')).astype(float))
    df['dist'], df['track_speed'] = gps_speed(df['lon'],df['lat'],df['DT'])
    df['track_direction'] = angles(df['lon'],df['lat'])
    return df.dropna().reset_index()[['DT','lon','lat','dt','dist','track_speed','track_direction']]

def A1inv(x):
    """Copy of A1inv function from circular package in R. Inverse function of the ratio of the first and zeroth order Bessel functions of the first kind. This function is used to compute the maximum likelihood estimate of the concentration parameter of a von Mises distribution.

    Args:

        x:  Value between 0 and 1

    Returns:
        Value of x such that A1inv(x) = k and A1(k) = x where A1 is the ratio between the first and zeroth order Bessel function of the first kind
    """
    if ((0 <= x) & (x < 0.53)):
        return 2 * x + pow(x,3) + (5 * pow(x,5))/6
    else:
        if (x < 0.85):
            return -0.4 + 1.39 * x + 0.43/(1-x)
        else:
            return 1/(pow(x,3) - 4 * pow(x,2) + 3 * x)

def rangeGen(DF,cen,wws,cutV,cutT,cutLength):
    """Test window for 

    Args:

        DF:         dataframe of tag recordings as output by prePare()
        cen:        index of the center value to be tested
        wws:        window width in seconds
        cutV:       minimum speed to ensure flight behaviour
        cutT:       maximum time difference between consecutive GPS points in seconds
        cutLength:  minimum number of samples allowed for model window

    Returns:
        If window conforms to model requirements (number of samples, flight speed, etc.), the range of indeces of the appropriate window and the center index of the window are returned
    """
    st = DF.dt[:cen][::-1].cumsum().gt(wws).idxmax() + 1
    end = DF.dt[(cen+1):].cumsum().gt(wws).idxmax() + 1
    if sum((DF.track_speed[st:end] > cutV) & (DF.dt[st:end] < cutT) & (DF.track_direction[st:end] != 100)) > cutLength:
        return range(st,end),cen

def findWindows(DF,cutv = 4.1667,windowlength = 51):
    """Calculate appropriate data windows for application of wind model

    Args:  

        DF:             dataframe of tag recordings as output by prePare()
        cutv:           minimum speed to ensure flight in m/s. Defaults to 4.1667
        windowlength:   window length duration in minutes. Defaults to 51

    Returns:
        Two arrays, a range array of windows that meet the model requirements, and an int array of center points of each window
    """

    # start from minimum possible point
    fs = (1/np.abs(np.timedelta64(mode(np.diff(DF.DT)),'m')).astype(int)).astype(int) # in fixes per minute
    expSamp = round(51 * fs) # expected number of samples
    cutlength = round(45/51 * expSamp)
    error_of_sampling_interval = 5 * fs # give 5 seconds of leeway for samples to be found in
    cutt = ((60 * fs) + error_of_sampling_interval).astype(int)
    windwidthsec = ((windowlength * fs)/2) * 60 + (60 / (4*fs))
    # go through each possible center value
    centr = DF.dt.cumsum().gt(windwidthsec).idxmax() + 1
    windows,centers = zip(*filter(lambda item: item is not None,[rangeGen(DF,center,windwidthsec,cutv,cutt,cutlength) for center in range(centr,DF.dt[::-1].cumsum().gt(windwidthsec).idxmax())]))
    return windows,centers

def initPars(head,spd,hed,cv = 34.7/3.6):
    """Generate initial parameters for log-likelihood testing

    Args:
        head:   heading direction (this value will range between -3:3 in testing)
        spd:    array of window track speeds
        hed:    array of window track headings
        cv:     an assumed constant air speed (m/s). Defaults to 34.7/3.6 as per Shiomi et al. 2012

    Returns:
        Mean heading (float), an initial value for kappa, the von Mises concentration parameter, and a length 5 array of initial values for log-likelihood testing
    """

    inita = 0
    while (inita < 5):
        inita = np.abs(np.random.normal(12.5,5))
    meangd = np.arctan2(np.mean(np.sin(hed)),np.mean(np.cos(hed)))
    inithd = meangd + (head/3 * np.pi/2)
    initkappa = A1inv(np.mean(np.cos(hed - meangd)))
    initmux = initkappa * np.cos(inithd)
    initmuy = initkappa * np.sin(inithd)
    initwx = np.mean(spd) * np.cos(meangd) - cv * np.cos(inithd)
    initwy = np.mean(spd) * np.sin(meangd) - cv * np.sin(inithd)
    return meangd,initkappa,[inita,initmux,initmuy,initwx,initwy]

def windOptims(spd,hed,cv,pars):
    """Optimisation for log-likelihood function produced by Likelihoodww() using Nelder-Mead methods

    Args:
        spd:    array of window track speeds
        hed:    array of window track headings
        cv:     an assumed constant air speed (m/s). Defaults to 34.7/3.6 as per Shiomi et al. 2012
        pars:   initial parameters to produce log-likelihood function

    Returns:
        Length 5 array of optimal parameters given initial values.
    """
    return minimize(Likelihoodww(spd,hed,cv),pars,bounds=([0.01,None],[None,None],[None,None],[None,None],[None,None]),method="L-BFGS-B")

def yokoTate(optimAns,cv):
    """Calculates von Mises SD and mean product along with Weibull standard deviation

    Args:
        optimAns:   output of windOptims for calculation of von Mises and Weibull characteristics
        cv:         an assumed constant air speed (m/s). Defaults to 34.7/3.6 as per Shiomi et al. 2012

    Returns:
        Product of von Mises SD and mean and the Weibull standard deviation, both using optimised values
    """
    yoko = Von_Mises_sd(np.sqrt(optimAns[2]*optimAns[2] + optimAns[1]*optimAns[1])) * Weibull_mean(optimAns[0],cv/sp.special.gamma(1+1/optimAns[0]))
    tate = Weibull_sd(optimAns[0], cv/sp.special.gamma(1 + 1/optimAns[0]))
    return yoko, tate

def ensureOptimConv(optimAns,spd,hed,cv):
    """Repeat maximum likelihood estimation to ensure convergence

    Args:

        optimAns:   optimised results of maximum likelihood from windOptims()
        spd:        window speed values (m/s)
        hed:        window heading values (rad)
        cv:         an assumed constant air speed (m/s). Defaults to 34.7/3.6 as per Shiomi et al. 2012

    Returns:
        Converged maximum likelihood estimates, along with their outputs of yokoTate()
    """

    newAns = optimAns
    pars = [newAns.x[0],newAns.x[1],newAns.x[2],newAns.x[3],newAns.x[4]]
    # substitute first parameter if estimated < 0
    if newAns.x[0] >= 0.01:
        newAns = windOptims(spd,hed,cv,pars)
        return newAns, yokoTate(newAns.x,cv)
    else:
        return optimAns, yokoTate(optimAns.x,cv)

def windOptim(initHead,spd,head,cv):
    """Full run of MLE for wind optimisation for suitable data window

    Args:
        initHead:   initial value for heading
        spd:        window track speeds (m/s)
        head:       window directions (rad)
        cv:         an assumed constant air speed (m/s). Defaults to 34.7/3.6 as per Shiomi et al. 2012

    Returns:
        If MLE has converged and values pass the required tests, optimised values and resultant von Mises and Weibull characteristics
    """

    # run first optimisation
    _,_,par = initPars(initHead,spd,head,cv)
    answ = windOptims(spd,head,cv,par)
    # SD of heading vector perpendicular (yoko) and along (tate) mean direction
    yoko,tate=yokoTate(answ.x,cv)
    # repeat MLE to ensure convergence
    if 'tate' in locals():
        answ,[yoko,tate] = ensureOptimConv(answ,spd,head,cv)
        return answ,yoko,tate
    else:
        return np.nan,np.nan,np.nan

def headSpdDir(spd,hed,answ):
    """Calculate speed and direction of heading from x and y components

    Args:
        spd:    window track speeds (m/s)
        hed:    window track headings (rad)
        answ:   optimised values from maximum likelihood

    Returns:
        Speed (nr, m/s) and direction (nd, rad) of heading vector
    """

    nvx = spd * np.cos(hed) - answ.x[3]
    nvy = spd * np.sin(hed) - answ.x[4]
    nr = np.sqrt(nvx**2 + nvy**2)
    nd = np.arctan2(nvy,nvx)
    return nr,nd

def GOFtests(hed,nr,nd,answ,yoko,tate,cv):
    """Test model output from windOptim()

    Args:
        hed:    window track headings (rad)
        nr:     heading vector speed
        nd:     heading vector direction
        answ:   optimised MLE estimates
        yoko:   von Mises characteristics
        tate:   Weibull characteristics
        cv:     an assumed constant air speed (m/s). Defaults to 34.7/3.6 as per Shiomi et al. 2012

    Returns:
        Array of three booleans indicating passing the model's three conditions
    """

    # mean track direction
    meangd = np.arctan2(np.mean(np.sin(hed)),np.mean(np.cos(hed)))
    nrp = sp.stats.kstest(nr,partial(sp.stats.weibull_min.cdf,c=answ.x[0],scale=cv/sp.special.gamma(1+1/answ.x[0])))
    # direction of heading vector
    mu = math.atan2(answ.x[2],answ.x[1])
    kappa = np.sqrt(pow(answ.x[1],2) + answ.x[2]*answ.x[2])
    ndp = sp.stats.kstest(nd,'vonmises',args=[kappa,mu])
    # correlation test between direction and speed of heading vector
    cnrnd = sp.stats.pearsonr(nr,nd)
    # describe conditions required
    cond1 = (yoko/tate) > 1
    cond2 = (np.cos(meangd)*np.cos(mu) + np.sin(meangd)*np.sin(mu)) > 0
    cond3 = (nrp.pvalue > 0.05) * (ndp.pvalue > 0.05) * (cnrnd[1] > 0.05)
    return [cond1,cond2,cond3]
    
def stringify(vars):
    """Turn array into string array

    Args:

        vars:   array of variables to be turned into string

    Returns:

        String array of same length as vars

    Example:
        >>> test = np.arange(0,50,5)
        >>> stringOfTest = stringify(test)
    """

    return [str(x) for x in vars]

def maxLikeWind(r,d,cv):
    """Generate optimal MLE estimation and return estimates if all three conditions are met

    Args:

        r:  track speeds (m/s)
        d:  track headings (rad)
        cv: an assumed constant air speed (m/s). Defaults to 34.7/3.6 as per Shiomi et al. 2012

    Returns:
        Array of length 3 containing estimate parameters from MLE
    """

    hd_try = 3
    answ_best = np.nan
    max_like = np.nan
    for id_hd in range(-hd_try,hd_try):
    
        # perform optimisation routine
        answ,yoko,tate = windOptim(id_hd,r,d,cv)

        if hasattr(answ,'success'):
            if (not np.isnan(tate)) & (answ.status==0):
        
                # calculate speeds and headings
                nr,nd = headSpdDir(r,d,answ)
        
                # GOF tests
                if np.isnan(max_like) & (np.prod(GOFtests(d,nr,nd,answ,yoko,tate,cv)) == 1):
                    max_like = answ.fun
                    answ_best = answ

                if (answ.fun > max_like) & (np.prod(GOFtests(d,nr,nd,answ,yoko,tate,cv)) == 1):
                    max_like = answ.fun
                    answ_best = answ

    return answ_best

def windEstimation(file, cutv: float = 4.1667, cv = 34.7/3.6, windowLength: int = 51, rescaleTime: bool = True, isBp: bool = True):
    """        wind estimation from bird GPS track

    Args:
        filename:       location of AxyTrek data (BiP or X Manager formatted)
        cutv:           minimum ground speed in m/s (default: 4.1667) 
        cv:             mean air speed in m/s (default: 34.7 kph)
        windowLength:   window size for estimation in mins (default: 51)
        rescaleTime:    should original data be rescaled to 1 fix/min (default: True)
        isBp:           boolean, is data formatted by BiP (T) or X Manager (F)

    Returns:
        Pandas dataframe of wind estimation methods. Columns are: centerpoint datetime ('Time'), lat ('Lat'), lon ('Lon'), mean bird heading during window (rad, east is 0. 'MeanHead'), wind vector x ('X'), wind vector y ('Y')
    """

    # # write the header of the outfile
    # with open(outfile,'w',newline='') as f:
    #     writer = csv.writer(f, delimiter = ',')
    #     writer.writerow(['Time','Lat','Lon','MeanHead','X','Y'])

    # initiate value
    out = []
    dat = prePare(file, convertToMin = rescaleTime, isBip = isBp)

    # generate windows over which estimation method will be run
    try:
        windows,centers = findWindows(dat,cutv,windowLength)
    except:
        return

    # max likelihood calculations for wind estimation
    for win in range(len(windows)):
        
        r = dat.track_speed[windows[win]]
        d = dat.track_direction[windows[win]]
        answ_best = maxLikeWind(r,d,cv)

        if type(answ_best) != float:
            out.append([(dat.DT[centers[win]]),dat.lat[centers[win]],dat.lon[centers[win]],np.arctan2(answ_best.x[2],answ_best.x[1]),answ_best.x[3],answ_best.x[4]])

    if len(out) > 0:
        out = pd.DataFrame(out)
        out.columns = ['Time','Lat','Lon','MeanHead','X','Y']
        return out

def roundTo(x,rn):
    return np.format_float_positional(x, precision=rn, unique=False, fractional=False, trim='k')

def makeForGraphText(num,rn):
    if (num < 0.1) & (num != 0):
        afterE = round(math.log10(num) - (0.5 if math.log10(num)<0 else 0))
        precede = roundTo(num * 10**abs(afterE),3)
        return str(precede) + "e" + str(afterE)
    else:
        return roundTo(num,rn)

# taken from 'https://gist.github.com/kn1cht/89dc4f877a90ab3de4ddef84ad91124e', circular python package by kn1cht
def Circcorrcoef(x, y, deg=True, test=False):
    '''Circular correlation coefficient of two angle data(default to degree)
    Set `test=True` to perform a significance test.
    '''
    convert = np.pi / 180.0 if deg else 1
    sx = np.frompyfunc(np.sin, 1, 1)((x - Circmean(x, deg)) * convert)
    sy = np.frompyfunc(np.sin, 1, 1)((y - Circmean(y, deg)) * convert)
    r = (sx * sy).sum() / np.sqrt((sx ** 2).sum() * (sy ** 2).sum())

    if test:
        l20, l02, l22 = (sx ** 2).sum(),(sy ** 2).sum(), ((sx ** 2) * (sy ** 2)).sum()
        test_stat = r * np.sqrt(l20 * l02 / l22)
        p_value = 2 * (1 - sp.stats.norm.cdf(abs(test_stat)))
        return tuple(round(v, 7) for v in (r, test_stat, p_value))
    return round(r, 7)

def Circmean(angles, deg=True):
    '''Circular mean of angle data(default to degree)
    '''
    a = np.deg2rad(angles) if deg else np.array(angles)
    angles_complex = np.frompyfunc(cmath.exp, 1, 1)(a * 1j)
    mean = cmath.phase(angles_complex.sum()) % (2 * np.pi)
    return round(np.rad2deg(mean) if deg else mean, 7)

# for BiP batch
def prePare2(df, convertToMin: bool = True, tdiff = 1, units = 'min', isBip: bool = True):
    """Prepare AxyTrek GPS data as per required for Goto original method. DateTime, distance, track speed and direction is added

    Args:
        filename:       full path for file to be read in
        convertToMin:   boolean for whether data should be resampled to regular time intervals
        tdiff:          resample time interval. Defaults to 1
        units:          resample time interval unit. Defaults to 'min'
        isBip:          boolean stating whether data is in BiP format

    Returns:
        Pandas dataframe of read in data. Adds columns 'dt' (elapsed time from fix to previous time point in seconds), 'dist' (distance travelled from previous point in m), 'track_speed' (in m/sec), 'track_direction' (direction between consecutive GPS positions in rad)
    """

    if convertToMin:
        df = timeRescale(df, tdiff, units)
    
    if len(df) > 1: 
        df['dt'] = np.append(np.nan,(np.diff(df['DT']) / np.timedelta64(1,'s')).astype(float))
        df['dist'], df['track_speed'] = gps_speed(df['lon'],df['lat'],df['DT'])
        df['track_direction'] = angles(df['lon'],df['lat'])
        return df.dropna().reset_index()[['DT','lon','lat','dt','dist','track_speed','track_direction']]
    else:
        return df

# for BiP batch
def windEstimation2(df, cutv: float = 4.1667, cv = 34.7/3.6, windowLength: int = 51, rescaleTime: bool = True, isBp: bool = True):
    """        wind estimation from bird GPS track

    Args:
        filename:       location of AxyTrek data (BiP or X Manager formatted)
        cutv:           minimum ground speed in m/s (default: 4.1667) 
        cv:             mean air speed in m/s (default: 34.7 kph)
        windowLength:   window size for estimation in mins (default: 51)
        rescaleTime:    should original data be rescaled to 1 fix/min (default: True)
        isBp:           boolean, is data formatted by BiP (T) or X Manager (F)

    Returns:
        Pandas dataframe of wind estimation methods. Columns are: centerpoint datetime ('Time'), lat ('Lat'), lon ('Lon'), mean bird heading during window (rad, east is 0. 'MeanHead'), wind vector x ('X'), wind vector y ('Y')
    """

    # # write the header of the outfile
    # with open(outfile,'w',newline='') as f:
    #     writer = csv.writer(f, delimiter = ',')
    #     writer.writerow(['Time','Lat','Lon','MeanHead','X','Y'])

    # initiate value
    out = []
    dat = prePare2(df, convertToMin = rescaleTime, isBip = isBp)
    
    if len(dat) >1:
        # generate windows over which estimation method will be run
        windows,centers = findWindows(dat,cutv,windowLength)

        # max likelihood calculations for wind estimation
        for win in range(len(windows)):
            
            r = dat.track_speed[windows[win]]
            d = dat.track_direction[windows[win]]
            answ_best = maxLikeWind(r,d,cv)

            if type(answ_best) != float:
                out.append([(dat.DT[centers[win]]),dat.lat[centers[win]],dat.lon[centers[win]],np.arctan2(answ_best.x[2],answ_best.x[1]),answ_best.x[3],answ_best.x[4]])

        if len(out) > 0:
            out = pd.DataFrame(out)
            out.columns = ['Time','Lat','Lon','MeanHead','X','Y']
            return out


def main_func(df):
    df['time'] = pd.to_datetime(df['time'], format='ISO8601')
    df = df.drop_duplicates('time', keep='last')
    df = df.dropna()
    df = df.reset_index()
    df = df.drop(['index'], axis=1)
    df.rename(columns={'time': 'DT', 'latitude': 'lat', 'longitude': 'lon'}, inplace=True)
    df['DT'] = [dtFormat(str(x)) for x in df['DT']]
    df['DT'] = pd.to_datetime(df['DT'],format="%Y-%m-%d %H:%M:%S.%f") #一旦、+00:00を除くために文字列にする
    
    if len(df) >0:
        out = windEstimation2(df, isBp = True)
        out = out.rename(columns={"Time":"time", "Lat":"latitude", "Lon":"longitude"})
        
        #すでにUTCなので不要
        out['time'] = pd.to_datetime(out['time'], format='%Y-%m-%d %H:%M:%S.%f', errors='coerce').dt.tz_localize('UTC')

        return out
    else:
        return pd.DataFrame(columns=['Time','Lat','Lon','MeanHead','X','Y'])
        