# DVLRead

from pathlib import Path

import pandas as pd
import re

# parent folder for all DVL data
dvlFolders = "E:/My Drive/PhD/Data/2018Shearwater/DVL"

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

    # add time series
    dat['DT'] = pd.date_range(tagInfo[re.search(r"(\d{5}?).txt",path.name).group(1)][0], periods = len(dat), freq = '200ms')

    # select data within video range
    dat = dat[(dat.DT >= tagInfo[re.search(r"(\d{5}?).txt",path.name).group(1)][1]) & (dat.DT < tagInfo[re.search(r"(\d{5}?).txt",path.name).group(1)][1] + pd.Timedelta(hours=2))]

    # read in the video behaviour records
    behavClass = pd.read_csv('E:/My Drive/PhD/Data/2018Shearwater/DVL/DiveBehavClass.csv', sep = ',', usecols = [0,1,2,4], dtype = {'Tag' : str})

    # replace 'Dive' behaviour with 's' (surface) or 'd' (dive)
    behavClass.loc[behavClass['Behaviour'] == 'Dive','Behaviour'] = behavClass.loc[behavClass.Behaviour == 'Dive','ForageBeh']

    # remove superfluous column
    behavClass.drop('ForageBeh',axis = 1, inplace = True)

    # extract behaviours for the specific tag
    behavClass.drop(behavClass[behavClass.Tag != re.search(r"(\d{5}?).txt",path.name).group(1)].index, inplace = True)

    tagData[re.search(r"(\d{5}?).txt",path.name).group(1)] = {'Acc' : dat, 'Beh' : behavClass}

# returns dictionary for each tag, each contain two keys, 'Acc' for acceleration data, 'Beh' for video-derived behaviours. Data should be aligned.