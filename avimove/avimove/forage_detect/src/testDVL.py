# DVLRead

from pathlib import Path
# bring in functions from main_func
from main_func import *
from math import pi

import pandas as pd
import re

# parent folder for all DVL data
# dvlFolders = "E:/My Drive/PhD/Data/2018Shearwater/DVL"
dvlFolders = "/Users/aran/Documents/DVL"

# acc and vid starttimes
st17008 = [pd.to_datetime('31/08/2018 06:00:00', format = "%d/%m/%Y %H:%M:%S"),
           pd.to_datetime('31/08/2018 12:00:00', format = "%d/%m/%Y %H:%M:%S")]
st18012 = [pd.to_datetime( '31/08/2018 04:00:00', format = "%d/%m/%Y %H:%M:%S"),
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

# read in the video behaviour records
behavClass = pd.read_csv('/Users/aran/Library/CloudStorage/GoogleDrive-aran.garrod@googlemail.com/My Drive/Data/DVL/DiveBehavClass.csv', sep = ',', usecols = [0,1,2,4], dtype = {'Tag' : str})

# replace 'Dive' behaviour with 's' (surface) or 'd' (dive)
behavClass.loc[behavClass['Behaviour'] == 'Dive','Behaviour'] = behavClass.loc[behavClass.Behaviour == 'Dive','ForageBeh']
# remove superfluous column
behavClass.drop('ForageBeh',axis = 1, inplace = True)

# set time to datetime
behavClass.Time = pd.to_datetime('31/08/2018 ' + behavClass.Time, format = '%d/%m/%Y %H.%M.%S.%f')


# path=Path(dvlFolders+"/18012/acc-video-2018-18012.txt")

for path in Path(dvlFolders).rglob('*acc*.txt'):
    
    # read in acceleration data
    dat = pd.read_table(path.__str__(), skiprows = 7, sep = ',', usecols = [0,1,2])

    # remove whitespace from headers
    dat.rename(columns=lambda x: x.strip(), inplace = True)

    # add time series
    dat['DT'] = pd.date_range(tagInfo[re.search(r"(\d{5}?).txt",path.name).group(1)][0], periods = len(dat), freq = '50ms')

    # select data within video range
    dat = dat[(dat.DT >= tagInfo[re.search(r"(\d{5}?).txt",path.name).group(1)][1]) & (dat.DT < tagInfo[re.search(r"(\d{5}?).txt",path.name).group(1)][1] + pd.Timedelta(hours=2))]

    # calculate pitch (clip array to min max of -1 1, assumes acceleration measured in g)
    dat = pd.concat([dat.reset_index(), accFeatures(dat[['X','Y','Z']],['Y','Z','X'],1.5,3,20)], axis = 1)
    tagname=re.search(r"(\d{5}?).txt",path.name).group(1)
    print(f"read in tag {tagname}")

    # generate spectrogram and summed spectral energy difference
    f,s,Sxx = hammingSpect(dat.Z, fs = 20)
    rollSum = rollingSpecSum(Sxx, f, 3, 5, fs = 20)
    print(f"Performing roll sum")

    # extract behaviours for the specific tag
    bc = behavClass[behavClass.Tag == re.search(r"(\d{5}?).txt",path.name).group(1)]
    print(f"Popping behaviour classes")

    # remove minute containing AT behaviour
    if bc.Behaviour.eq('AT').any():
        print(f"AT present")
        behAT = np.any([(dat.DT >= (x - pd.Timedelta(30,'sec'))) & (dat.DT <= (x + pd.Timedelta(30,'sec'))) for x in bc.Time[bc.index[bc.Behaviour == 'AT']].round("s").values], axis = 0)

        rollSum[behAT[(2*20):-(2*20)+1]] = min(rollSum)

    tagData[re.search(r"(\d{5}?).txt",path.name).group(1)] = {'acc' : dat, 'beh' : bc}

# returns dictionary for each tag, each contain two keys, 'acc' for acceleration data, 'beh' for video-derived behaviours. Data should be aligned.

maxWithGap(rollSum,20,numPoints=10)

# fs=20
# window=60
# [np.arange(np.argmax(sigad) - round(fs*window/2),
#     np.argmax(sig) + round(fs*window/2))]

# out = []
# sigad = rollSum.copy()
# while len(out) != 5:
#     # create index range around highest value
#     out.append([np.arange(np.argmax(sigad) - round(fs*window/2),
#                 np.argmax(sigad) + round(fs*window/2))])
#     # reduce magnitude of this period and 5 minutes surrounding
#     sigad[np.arange(np.argmax(sigad) - round(fs*window/2) - (fs*60*5),
#                 np.argmax(sigad) + round(fs*window/2)) + (fs*60*5)] = np.min(sigad)
# out

# np.argmax(sigad) - round(fs*window/2)

# tagData['17008']['beh']

# dat.iloc[:,[0,4]].diff(axis=1).iloc[:,1].values

# range(2).apply(lambda x: dat.iloc[:,x] - dat.iloc[:,x + 4])

# dat.iloc[:,range(2)] - dat.iloc[:,range(4,7)]


# dat.iloc[:,0]


# # for tag in tagData:
# #     tag['acc'].loc[:,['SY','SX','SZ']]

# tagData['17008']['acc'].loc['SX','SY','SZ'] = 

# tagData['17008']['acc'].assign([lowEquiFilt(x,1,1.5,25) for x in zip([tagData['17008']['acc'].X,tagData['17008']['acc'].Y,tagData['17008']['acc'].Z])])



# df[cols]=df[cols].apply(lambda x: pd.to_datetime(x ,errors='coerce')).dt.strftime('%Y-%b-%d')


# for tag in tagData.keys():
#     static = dict(zip(['SY','SX','SZ'],))