import re

import numpy as np
import pandas as pd
import geopy.distance as distance
from scipy import stats

def readAxyGPS(filename, delim = "\t", cols = None, datetimeFormat = "%d/%m/%Y %H:%M:%S"): 
    """
    Read in AxyTrek GPS data (txt files) as output by X Manager

    Args:

        filename:   path to AxyTrek txt file
        cols:       string to denote if acceleration and GPS ('acc') or GPS only ('gps') should be read in
        colnames:   list of string column names to be assigned. Must be of same length as cols. Defaults to ['Date', 'Time', 'lat', 'lon']

    Returns:
        Pandas dataframe of columns `colnames`. A formatted DateTime column (named DT) is generated.
    """

    colnames = ['Date', 'Time', 'lat', 'lon']


    df = pd.read_csv(filename, sep = delim, usecols = cols,
    names = colnames)
    df.DT = [dtFormat(x) for x in df.DT] # ensure correct datetime formats
    df['DT'] = pd.to_datetime(df['Date'] + " " + df['Time'], format = datetimeFormat)
    return df

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

def distSpeed(lat,lon,DT):

    dist,speed = gps_speed(lat,lon,DT)
    while np.nanmax(speed) > threshold:

        # remove erroneous GPS values
        lon[speed > threshold] = np.nan
        lat[speed > threshold] = np.nan

        # recalculate speed
        dist,speed = gps_speed(lat,lon,DT)
    
    return dist,speed

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

dat = readBIP("https://bipsharedata.s3.ap-northeast-1.amazonaws.com/analysis/cal_wind2/Axy/3a66a690-5d53-4bed-a4a7-4d25a1569c4d/3a66a690-5d53-4bed-a4a7-4d25a1569c4d_std.csv",cols='acc')

tst = dat.iloc[0:1610000,].reset_index()
tstGPS = tst.dropna()

(110000-20000)/(25*60)

# remove positions near capture site
FkOshima = [39.400,141.998] # capture site
test = np.array(gps_distanceSingle(dat.lon.dropna(),dat.lat.dropna(),FkOshima[0],FkOshima[1]))
sum(test<1500)

test < 15
test[-1]



np.where(test < 1500)[0][-1]
import matplotlib.pyplot as plt

plt.scatter(dat.lat.dropna(),dat.lon.dropna())
plt.scatter(dat.lat.dropna()[test < 5],dat.lon.dropna()[test < 5],color='red')
plt.show()
# find GPS sampling rate
np.timedelta64(stats.mode(np.diff(tst.DT.dropna()))[0][0],'ms')

sum(test < 1.5)

rem = []
for x in np.where(test < 1.5)[0]:
    rem.extend(np.arange(np.max([0,x-int((60/gpsFS * fs)/2)]),np.min([len(dat),x+int((60/gpsFS * fs)/2)])))
rem = np.array(rem)

dat.drop(rem)

b = np.where(test < 1.5)[0][0]
np.arange(np.max([0,b-int((60/gpsFS * fs)/2)]),np.min([len(dat),b+int((60/gpsFS * fs)/2)]))

tst.dropna()

# acceleration sample rate (Hz)
fs = int(np.timedelta64(1,'s') / stats.mode(np.diff(dat.DT))[0][0])

# GPS sample rate (fixes per minute)
gpsFS = int(np.timedelta64(1,'m') / stats.mode(np.diff(dat.dropna()['DT']))[0][0])


test[test != np.nan]

test[0] == nan

len(np.where(test < 1500)[0])

len(tst.dropna())




1/(np.timedelta64(1,'s')/stats.mode(np.diff(dat.dropna().DT))[0][0])

gpsMode = stats.mode()
lats = 
lons = dat.lon.dropna()


for x in np.where(test < 1500)[0]:
    tst.drop(range(np.max([0,x-12]),np.min([x+12,len(tst)-1])), inplace = True)

tst.drop(range(90000,95000), axis=0)



rem = rem[rem[0] < len(tst)]

rem > list(len(tst))

type(rem)

type(len(tst))


rem < len(tst)

tst.drop(rem,axis=1)

np.max(rem)

[distance.distance((x,y),FkOshima).m for x,y in zip(tstGPS.lat,tstGPS.lon)]