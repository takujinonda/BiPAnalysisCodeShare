# Script to assess accuracy of behaviour detection method on DVL tag data

import platform 
import pandas as pd

from birdTag import birdTag
from pathlib import Path

import glob

# parent folder for all DVL data
if platform.system() == 'Darwin':
    dvlFolders = "/Users/aran/Documents/DVL/"
else:
    dvlFolders = "E:/My Drive/PhD/Data/2018Shearwater/DVL/"

# read in the video behaviour records
if platform.system() == 'Darwin':
    behav_class_loc = '/Users/aran/Library/CloudStorage/GoogleDrive-aran.garrod@googlemail.com/My Drive/Data/DVL/DiveBehavClass.csv'
else:
    behav_class_loc = 'E:/My Drive/PhD/Data/2018Shearwater/DVL/DiveBehavClass.csv'

accFiles = Path(dvlFolders).rglob('*/acc*.txt')

tags = [
    {'name':'17008','accSt':'31/08/2018 06:00:00','vidSt':'31/08/2018 12:00:00'},
    {'name':'18012','accSt':'31/08/2018 04:00:00','vidSt':'31/08/2018 10:00:00'},
    {'name':'18014','accSt':'31/08/2018 06:00:00','vidSt':'31/08/2018 12:00:00'},
    {'name':'18017','accSt':'31/08/2018 04:00:00','vidSt':'31/08/2018 11:00:00'},
    {'name':'18018','accSt':'31/08/2018 04:00:00','vidSt':'31/08/2018 12:00:00'}
]

dvls = []
for ind,tag in enumerate(tags):
    dvls.append({'name':tag['name'],'data':birdTag(filepath=glob.glob(dvlFolders+tag['name']+'/*acc*.txt')[0],type='dvl',tagname=tag['name'],accfs=20,long_acc_name=['Y','Z','X'],accStart=tag['accSt'],vidStart=tag['vidSt'])})
    dvls[ind]['data'].readin()
    dvls[ind]['data'].readBeh(behavPath=behav_class_loc)

dvls[0]['data'].flight_est(10,True)

dvls[0]['data'].readBeh(behavPath=behav_class_loc)