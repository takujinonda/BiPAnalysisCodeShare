from pathlib import Path
# bring in functions from main_func
from forage_detect.main_func import *
from math import pi

import pandas as pd
import re
import argparse

class birdTag:
    """Class for avian biologging tags. Current support includes AxyTrek and DVL with a look toward incorporating NinjaScan data.
    Requires file path(s) and string descriptor of tag ('Axy','DVL').
    Methods include reading of data from raw CSV/txt or BiP database. Data loads can include GPS or not (where present).
    Vectorised haversine distance and speed calculation from GPS.
    Dynamic/static acceleration calculation via equiripple lowpass filter.
    Pitch calculation.
    Data removal from proximity to location.
    """
    def __init__(self,filepath,type,tagname,fs):
        self.filepath = filepath
        self.type = type
        self.tagname = tagname
        self.fs = fs

    def readAxyGPS(self, delim = "\t", cols = None, datetimeFormat = "%d/%m/%Y %H:%M:%S"): 
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


        df = pd.read_csv(self.filename, sep = delim, usecols = cols,
        names = colnames)
        df.DT = [self.dtFormat(x) for x in df.DT] # ensure correct datetime formats
        df['DT'] = pd.to_datetime(df['Date'] + " " + df['Time'], format = datetimeFormat)
        df.rename(columns={'Date':'date','Time':'time','DT':'dt'},inplace=True) # standardise columns
        self.gps = df
    
    def readInAxyCSV(self,sep='\t',cols=None):
        if cols is None:
            cols = ['TagID','Timestamp','X','Y','Z','location-lat','location-lon']
        df = pd.read_csv(self.filename,usecols=cols,sep=sep)
        df.loc[:,'Timestamp'] = pd.to_datetime(df.loc[:,'Timestamp'],format="%d/%m/%Y %H:%M:%S.%f")
        df.rename(columns={'TagID':'tagid','Timestamp':'dt','location-lat':'lat','location-lon':'lon'},inplace=True)
        self.acc = df[['tagid','dt','X','Y','Z']]
        self.gps = df[['tagid','dt','lat','lon']].dropna()
    
    def readAxyBIP(self,cols=None):
        """
        Read in BiP-formatted AxyTrek data

        Args:

            filename:   string of full path to file
            cols:       string value indicating if number of columns to read in can be reduced. For GPS and acceleration, use 'acc', for GPS only, use 'gps'. Defaults to None.

        Returns:
            Pandas dataframe of BIP system formatted Axytrek data. Columns returned depend on `cols` argument. 'acc' returns formatted datetime (DT), latitude (lat), longitude (lon), longitudinal acceleration (X), lateral acceleration (Y), and dorsoventral acceleration (Z). 'gps' returns datetime, latitude, and longitude. If no cols argument given, all data is read in, including pressure (pressure), temperature (temp), height about sea level (altitude), and ground speed (spd).
        """

        accCols = ['time','latitude','longitude','acceleration_longitudinal','acceleration_lateral','acceleration_dorso_ventral']
        gpsCols = ['time','latitude','longitude']

        # read in based on requested columns
        if cols.lower() == 'gps':
            df = pd.read_csv(self.filename, sep = ",", header = 0, usecols = gpsCols)
        elif cols.lower() == 'acc':
            df = pd.read_csv(self.filename, sep = ",", header = 0, usecols = accCols)
        else:
            df = pd.read_csv(self.filename, sep = ",", header = 0)
            
        # df = pd.read_csv(filename, sep = ",", header = 0, usecols = cols)
        # rename columns for later use
        df.rename(columns = {'time':'dt','latitude':'lat','longitude':'lon','acceleration_longitudinal':'X','acceleration_lateral':'Y','acceleration_dorso_ventral':'Z','pressure':'pressure','temperature':'temp','height_above_mean_sea_level':'altitude','ground_speed':'spd'}, inplace = True)
        
        df.DT = [self.dtFormat(x) for x in df.DT] # ensure correct datetime formats
        df['DT'] = pd.to_datetime(df['DT'], format = "%Y-%m-%d %H:%M:%S.%f")
        self.acc = df[['dt','X','Y','Z']]
        self.gps = df[['dt','lat','lon']].dropna()
    
    def readDVL(self,accStart,vidStart,removePriorVid=True):
        """
        Read Little Leonardo DVL data logger (400M) data. File locations and start-times of acceleration and video recording (dd/mm/yyyy HH:MM:SS) required.
        """
        # read in acceleration data
        self.acc = pd.read_table(self.filepath, skiprows=7, sep=',', usecols=[0,1,2])
        # remove header whitespace
        self.acc.rename(columns=lambda x: x.strip(), inplace=True)
        self.dat['DT'] = pd.date_range(accStart, periods=len(self.dat), freq=f"{1/self.fs * 10**3}ms")# add time series
        if removePriorVid:
            # select data within video range
            self.dat = self.dat[(self.dat.DT >= vidStart) & (self.dat.DT < (vidStart + pd.Timedelta(hours=2)))]
        
    
    def gps_distanceSingle(self, latVal, lonVal):                                                                       
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

        return [distance.distance((x,y),(latVal,lonVal)).km for x,y in zip(self.gps.lat,self.gps.lon)]
    
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

    def gps_speed(self):                                                                       
    # taken from https://www.tjansson.dk/2021/03/vectorized-gps-distance-speed-calculation-for-pandas/, thanks to Thomas Jansson for use of this function
        """                                                                                                                 
        Calculates the instantaneous speed from the GPS positions and timestamps. The distances between the points          
        are calculated using a vectorized haversine calculation the great circle distance between two arrays of points on   
        the earth (specified in decimal degrees). All args must be of equal length.                                         
    
        Returns:                                                                                                            
            Speed is returned an array in m/s.                                                                             
    
        Example:                                                                                                            
            >>> df['gpsSpeed'] = gps_speed(df.longitude, df.latitude, df.recordedAt).
        """
    
        lon1 = self.gps.lon[:-1].reset_index(drop = True)                                                                                       
        lat1 = self.gps.lat[:-1].reset_index(drop = True)                                                                                 
        lon2 = self.gps.lon[1:].reset_index(drop = True)                                                                                       
        lat2 = self.gps.lat[1:].reset_index(drop = True)                                                                                 
    
        # Vectorized haversine calculation                                                                                  
        lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])                                                  
        a = np.sin((lat2 - lat1) / 2.0)**2 + (np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2.0)**2)                 
        dist = np.array(6371 * 2 * np.arcsin(np.sqrt(a)) * 1000)
        time_array = (self.gps.dt.diff().dt.seconds ).values[1:]                                                       
    
        # Calculate the speed                                                                                               
        time_array[time_array == 0] = np.nan  # To avoid division by zero                                                   
        speed = np.array(dist / time_array)

        # Make the arrays as long as the input arrays                                                                        
        self.speed = np.insert(speed, 0, np.nan, axis=0)
        self.dist = np.insert(dist, 0, np.nan, axis=0)                                                               

    def distSpeed(self,threshold=None):
        """
        Calculates distances (metres) and speed (m/s) of GPS positions and datetime information. Uses a speed threshold (default None) to define erroneous GPS positions. GPS and datetime lengths must agree.

        Args:

            lat:    array of latitude points in decimal degrees.
            lon:    array of longitude points in decimal degrees.
            DT:     datetime array in Pandas datetime format.

        Returns:
            Float arrays of distance and speed values.
        """

        dist,speed = gps_speed(self.gps.lat,self.gps.lon,self.gps.dt)
        if threshold is not None:

            while np.nanmax(speed) > threshold:

                # remove erroneous GPS values
                self.gps.lon[speed > threshold] = np.nan
                self.gps.lat[speed > threshold] = np.nan

                # recalculate speed
                dist,speed = gps_speed(self.gps.lat,self.gps.lon,self.gps.dt)
        
        self.speed = speed
        self.dist = dist

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

    def removeNear(self, data, lat, lon, distThreshold = 1.5):
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
        fs = int(np.timedelta64(1,'s') / stats.mode(np.diff(self.acc.dt),keepdims=False)[0])

        # GPS sample rate (fixes per minute)
        gpsFS = 60/int(np.timedelta64(1,'m') / stats.mode(np.diff(data.gps.dt),keepdims=False)[0])

        # find distances from capture site
        capSiteDist = np.array(gps_distanceSingle(lat,lon))

        inds = data.gps.lat.index[np.where(capSiteDist < distThreshold)[0]] # find indeces to be removed from original dataset

        # if any difference in GPS and acc sampling rates means we must extend GPS indeces
        samplerDiff = int(np.floor(fs * gpsFS / 2)) # number of fs samples either side of GPS position

        # create ranges to be removed
        rem = []
        for x in inds:
            rem.extend(np.arange(x - samplerDiff, x + samplerDiff))
        rem = np.array(rem)

        self.acc.drop(rem,inplace=True)
        self.gps.drop(inds,inplace=True)

    def lowEquiFilt(self,sig, passband, stopband):
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
        eqFil=signal.remez(101,[0,passband,stopband,fs*.5],[1,0],fs=self.fs)
        # return 'static' signals for each provided signal
        static = signal.filtfilt(b=eqFil,a=1,x=sig)
        # return 'dynamic' acceleration signal
        dynamic = np.array(sig - static)
        
        return static, dynamic

    # create static and dynamic (1 pass, 1.5 stop), pitch, ODBA
    def accFeatures(self, passb, stopb):
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
        self.acc.sZ, self.acc.dZ = lowEquiFilt(self.acc.Z, passb, stopb, fs)
        self.acc.sX, self.acc.dX = lowEquiFilt(self.acc.X, passb, stopb, fs)
        self.acc.sY, self.acc.dY = lowEquiFilt(self.acc.Y, passb, stopb, fs)

        self.acc.pitch = np.arcsin(np.clip(sLong,-1,1)) * 180/np.pi
        self.acc.ODBA = sum([np.abs(self.acc.dZ),np.abs(self.acc.dX),np.abs(self.acc.dY)])

        self.acc.pitmn = self.acc.pitch.rolling(10*self.fs, closed = "both", min_periods = 1).mean()
        self.acc.ODmn = self.acc.ODBA.rolling(10*self.fs, closed = "both", min_periods = 1).mean()

    def hammingSpect(self,sig):
        """
        Generate spectrogram data of sig using a Hamming window of 4 seconds with 85% overlap.
        
        Args:
            sig:    signal to generate spectrogram.
            fs:     sig sampling frequency (Hz). Defaults to 25.
            
        """
        # generate spectrogram (4 second window, overlap of 85%, hamming window)

        #set window size of 4 seconds
        winsize = self.fs*4
        #set overlap between windows (will determine the temporal resolution
        numoverlap = np.floor(.9875*winsize).astype(int); #(85%)
        win = signal.windows.hamming(winsize);

        f, t, Sxx = signal.spectrogram(sig,self.fs,window=win,noverlap=numoverlap)

        return f, t, Sxx

    def rollingSpecSum(self,spec,f,minFreq,maxFreq,dur=60,inclusive=False):
        """
        Create a rolling sum of the difference between spectrogram intensities from minFreq to maxFreq and intensities above maxFreq, , i.e. summed values between 3 and 5 Hz - summed values above 5 Hz.
        
        Args:
            spec:       spectrogram object.
            f:          array of sample frequencies.
            minFreq:    minimum frequency for sum (Hz).
            maxFreq:    maximum frequency for sum (Hz).
            dur:        rolling sum duration (s).
            inclusive:  should boundaries be included. Defaults to False.

        Returns:
            Summed frequency intensities across spec object.
        """
        if inclusive:
            spectDiff = pd.Series(np.sum(spec[(f >= minFreq) & (f <= maxFreq),:],axis=0) - np.sum(spec[f > maxFreq,:],axis=0))
        else:
            spectDiff = pd.Series(np.sum(spec[(f > minFreq) & (f < maxFreq),:],axis=0) - np.sum(spec[f > maxFreq,:],axis=0))
        return spectDiff.rolling(dur*self.fs, closed = "both", min_periods = 1).sum()

    def maxWithGap(self,sig,window=60,minGap=5,numPoints = 20):
        """
        Calculate the highest values across time-series signal `sig` with a minimum time gap between
        
        Args:
            sig:        signal whose largest values are to be found.
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
            out.append([np.arange(np.argmax(sigad) - round(self.fs*window/2),
                    np.argmax(sig) + round(self.fs*window/2))])
            # reduce magnitude of this period and 5 minutes surrounding
            sigad[np.arange(np.argmax(sigad) - round(self.fs*window/2) - (self.fs*60*5),
                    np.argmax(sig) + round(self.fs*window/2)) + (self.fs*60*5)] = np.min(sigad)
        return out
    
    def DVLbehaviours(self,fileloc):
        self.behavClass = pd.read_csv(fileloc, sep=',', usecols=[0,1,2,4], dtype={'Tag': str})
        # replace 'Dive' behaviour with 's' (surface) or 'd' (dive)
        self.behavClass.loc[self.behavClass['Behaviour'] == 'Dive','Behaviour'] = self.behavClass.loc[self.behavClass.Behaviour == 'Dive','ForageBeh']

        # remove superfluous column
        self.behavClass.drop('ForageBeh',axis = 1, inplace = True)

        # extract behaviours for the specific tag
        self.behavClass.drop(self.behavClass[self.behavClass.Tag != self.tagname].index, inplace = True)

    def ATremove(self,rollSum):
        if not hasattr(self,'behavClass'):
            raise Exception('No behaviour class present (DVLbehaviours())')
        
        if self.behavClass.Behaviour.eq('AT').any():
            behAT = np.any([(self.acc.dt >= (x - pd.Timedelta(30,'sec'))) & (self.acc.dt <= (x + pd.Timedelta(30,'sec'))) for x in self.behavClass.Time[self.behavClass.index[self.behavClass.Behaviour == 'AT']].round("s").values], axis = 0)

        rollSum[behAT[(2*20):-(2*20)+1]] = min(rollSum)
