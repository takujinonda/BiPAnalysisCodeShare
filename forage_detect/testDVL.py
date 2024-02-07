# DVLRead

from pathlib import Path
# bring in functions from main_func
from forage_detect.main_func import *
from math import pi

import pandas as pd
import re

# parent folder for all DVL data
dvlFolders = "F:/My Drive/PhD/Data/2018Shearwater/DVL"

# acc and vid starttimes
st17008 = [pd.to_datetime('31/08/2018 06:00:00', format = "%d/%m/%Y %H:%M:%S"),
           pd.to_datetime('31/08/2018 12:00:00', format = "%d/%m/%Y %H:%M:%S")]
st18012 = [pd.to_datetime('31/08/2018 04:00:00', format = "%d/%m/%Y %H:%M:%S"),
            pd.to_datetime('31/08/2018 10:00:00', format = "%d/%m/%Y %H:%M:%S")]
st18014 = [pd.to_datetime('31/08/2018 06:00:00', format = "%d/%m/%Y %H:%M:%S"),
            pd.to_datetime('31/08/2018 12:00:00', format = "%d/%m/%Y %H:%M:%S")]
st18017 = [pd.to_datetime('31/08/2018 04:00:00', format = "%d/%m/%Y %H:%M:%S"),
            pd.to_datetime('31/08/2018 11:00:00', format = "%d/%m/%Y %H:%M:%S")]
st18018 = [pd.to_datetime('31/08/2018 04:00:00', format = "%d/%m/%Y %H:%M:%S"),
            pd.to_datetime('31/08/2018 12:00:00', format = "%d/%m/%Y %H:%M:%S")]
starttimes = [st17008,st18012,st18014,st18017,st18018]

tagnames = ['17008','18012','18014','18017','18018']

# dictionaries storing start times (acceleration [0], video [1])
tagInfo = dict(zip(tagnames,starttimes))

# initialise dictionary
tagData = {}

for path in Path(dvlFolders).rglob('*acc*.txt'):

    # read in acceleration data
    dat = pd.read_table(path.__str__(), skiprows = 7, sep = ',', usecols = [0,1,2])

    # remove whitespace from headers
    dat.rename(columns=lambda x: x.strip(), inplace = True)

    # add time series
    dat['DT'] = pd.date_range(tagInfo[re.search(r"(\d{5}?).txt",path.name).group(1)][0], periods = len(dat), freq = '200ms')

    # select data within video range
    dat = dat[(dat.DT >= tagInfo[re.search(r"(\d{5}?).txt",path.name).group(1)][1]) & (dat.DT < tagInfo[re.search(r"(\d{5}?).txt",path.name).group(1)][1] + pd.Timedelta(hours=2))]

    # add static and dynamic accelerations
    dat.loc[:,['sX','sY','sZ','dX','dY','dZ']] = np.concatenate(dat.loc[:,['X','Y','Z']].apply(lambda x: lowEquiFilt(x,1,1.5,20)).values, axis = 0)

    # calculate pitch (clip array to min max of -1 1, assumes acceleration measured in g)
    dat.loc[:,'pitch'] = np.arcsin(np.clip(dat.sY,-1,1)) * 180/pi

    # generate spectrogram and summed spectral energy difference
    f,s,Sxx = hammingSpect(dat.Z, fs = 20)
    rollSum = rollingSpecSum(Sxx, f, 3, 5, fs = 20)

    # read in the video behaviour records
    behavClass = pd.read_csv('F:/My Drive/PhD/Data/2018Shearwater/DVL/DiveBehavClass.csv', sep = ',', usecols = [0,1,2,4], dtype = {'Tag' : str})

    # replace 'Dive' behaviour with 's' (surface) or 'd' (dive)
    behavClass.loc[behavClass['Behaviour'] == 'Dive','Behaviour'] = behavClass.loc[behavClass.Behaviour == 'Dive','ForageBeh']

    # remove superfluous column
    behavClass.drop('ForageBeh',axis = 1, inplace = True)

    # extract behaviours for the specific tag
    behavClass.drop(behavClass[behavClass.Tag != re.search(r"(\d{5}?).txt",path.name).group(1)].index, inplace = True)

    # set time to datetime
    behavClass.Time = pd.to_datetime('31/08/2018 ' + behavClass.Time, format = '%d/%m/%Y %H.%M.%S.%f')

    # remove periods containing AT behaviour
    if behavClass.Behaviour.eq('AT').any():
        for b in behavClass.Time[behavClass.Behaviour == 'AT'].round('s').values:
            

mask = (dat.DT > find[0] - pd.Timedelta(30,'sec')) & (dat.DT < find[0] + pd.Timedelta(30,'sec'))

dat.DT[mask]
    
find = behavClass.Time[behavClass.index[behavClass.Behaviour == 's']].round("s").values


    tagData[re.search(r"(\d{5}?).txt",path.name).group(1)] = {'acc' : dat, 'beh' : behavClass}

# returns dictionary for each tag, each contain two keys, 'acc' for acceleration data, 'beh' for video-derived behaviours. Data should be aligned.

tagData['17008']
dat.loc[:,[['X','sX'],['Y','sY']]]


dat.iloc[:,[0,4]].diff(axis=1).iloc[:,1].values

range(2).apply(lambda x: dat.iloc[:,x] - dat.iloc[:,x + 4])

dat.iloc[:,range(2)] - dat.iloc[:,range(4,7)]


dat.iloc[:,0]


# for tag in tagData:
#     tag['acc'].loc[:,['SY','SX','SZ']]

tagData['17008']['acc'].loc['SX','SY','SZ'] = 

tagData['17008']['acc'].assign([lowEquiFilt(x,1,1.5,25) for x in zip([tagData['17008']['acc'].X,tagData['17008']['acc'].Y,tagData['17008']['acc'].Z])])



df[cols]=df[cols].apply(lambda x: pd.to_datetime(x ,errors='coerce')).dt.strftime('%Y-%b-%d')


for tag in tagData.keys():
    static = dict(zip(['SY','SX','SZ'],))