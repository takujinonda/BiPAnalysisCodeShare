# DVLRead

from pathlib import Path
# bring in functions from main_func
# from forage_detect.main_func import *
from math import pi

import pandas as pd
import re
import argparse

class DVL:
  """
  Class for Little Leonardo DVL data logger data. File locations and start-times of acceleration and video recording (dd/mm/yyyy HH:MM:SS)
  """
  def __init__(self, filepath, accStart, vidStart):
    try:
      self.file = filepath
      self.accStart = pd.to_datetime(accStart, format='%d/%m/%Y %H:%M:%S')
      self.vidStart = pd.to_datetime(vidStart, format='%d/%m/%Y %H:%M:%S')
    except ValueError as e:
      print(e)

  def download(self,removePriorVid=True):
    """
    Download data from filepath supplied 
    """
    # read in acceleration data
    self.dat = pd.read_table(self.filepath, skiprows=7, sep=',', usecols=[0,1,2])
    # remove header whitespace
    self.dat.rename(columns=lambda x: x.strip(), inplace=True)
    self.dat['DT'] = pd.date_range(self.accStart, periods=len(self.dat), freq='50ms')# add time series
    # select data within video range
    self.dat = self.dat[(self.dat.DT >= self.vidStart) & (self.dat.DT < (self.vidStart + pd.Timedelta(hours=2)))]
    return

  def pitchcal(self):
    """
    Generate pitch values (clips acceleration values to [-1,1])
    """
    self.pitch = pd.concat([self.dat.reset_index(), (self.dat[['X','Y','Z']],['Y','Z','X'],1.5,3,20)], axis=1, inplace=True)

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

    # remove whitespace from headers
    dat.rename(columns=lambda x: x.strip(), inplace = True)

    # add time series
    dat['DT'] = pd.date_range(tagInfo[re.search(r"(\d{5}?).txt",path.name).group(1)][0], periods = len(dat), freq = '50ms')

    # select data within video range
    dat = dat[(dat.DT >= tagInfo[re.search(r"(\d{5}?).txt",path.name).group(1)][1]) & (dat.DT < tagInfo[re.search(r"(\d{5}?).txt",path.name).group(1)][1] + pd.Timedelta(hours=2))]

    # calculate pitch (clip array to min max of -1 1, assumes acceleration measured in g)
    dat = pd.concat([dat.reset_index(), accFeatures(dat[['X','Y','Z']],['Y','Z','X'],1.5,3,20)], axis = 1, inplace=True)
    
    # generate spectrogram and summed spectral energy difference
    f,s,Sxx = hammingSpect(dat.Z, fs = 20)
    rollSum = rollingSpecSum(Sxx, f, 3, 5, fs = 20)

    # read in the video behaviour records
    behavClass = pd.read_csv('E:/My Drive/PhD/Data/2018Shearwater/DVL/DiveBehavClass.csv', sep = ',', usecols = [0,1,2,4], dtype = {'Tag' : str})

    # replace 'Dive' behaviour with 's' (surface) or 'd' (dive)
    behavClass.loc[behavClass['Behaviour'] == 'Dive','Behaviour'] = behavClass.loc[behavClass.Behaviour == 'Dive','ForageBeh']

    # remove superfluous column
    behavClass.drop('ForageBeh',axis = 1, inplace = True)

    # extract behaviours for the specific tag
    behavClass.drop(behavClass[behavClass.Tag != re.search(r"(\d{5}?).txt",path.name).group(1)].index, inplace = True)

    # set time to datetime
    behavClass.Time = pd.to_datetime('31/08/2018 ' + behavClass.Time, format = '%d/%m/%Y %H.%M.%S.%f')

    # remove minute containing AT behaviour
    if behavClass.Behaviour.eq('AT').any():
        behAT = np.any([(dat.DT >= (x - pd.Timedelta(30,'sec'))) & (dat.DT <= (x + pd.Timedelta(30,'sec'))) for x in behavClass.Time[behavClass.index[behavClass.Behaviour == 'AT']].round("s").values], axis = 0)

        rollSum[behAT[(2*20):-(2*20)+1]] = min(rollSum)

    tagData[re.search(r"(\d{5}?).txt",path.name).group(1)] = {'acc' : dat, 'beh' : behavClass}

# returns dictionary for each tag, each contain two keys, 'acc' for acceleration data, 'beh' for video-derived behaviours. Data should be aligned.
