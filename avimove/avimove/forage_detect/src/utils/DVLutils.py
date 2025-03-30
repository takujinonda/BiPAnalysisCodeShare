import pandas as pd
import numpy as np

from pathlib import Path
from .loadIn import readDVL

def readBeh(behaviour_file, name):
    out = pd.read_csv(behaviour_file, sep = ',', usecols = [0,1,2,4], dtype = {'Tag': str})
    # contain only the data of tag of interest
    out = out.loc[out.Tag == name].reset_index(drop=True)
    # replace 'Dive' behaviour with 's' (surface) or 'd' (dive)
    out.loc[out.Behaviour == 'Dive','Behaviour'] = out.loc[out.Behaviour == 'Dive','ForageBeh']
    # remove superfluous column
    out.drop('ForageBeh',axis=1,inplace=True)
    # change 'IT' behaviour (intermediate takeoff) to flight
    out.loc[out.Behaviour == 'IT','Behaviour'] = "FL"

    # set time to datetime
    out.Time = pd.to_datetime('31/08/2018 ' + out.Time, format = '%d/%m/%Y %H.%M.%S.%f')

    return out

def readAllDVL(location,names,stTimes,vidTimes):
    test = {}
    for i,path in Path(location).rglob('*/acc*.txt'):
        test.update({names[i]:{'acc':readDVL(path,stTimes[i],fs=20,vidStart=vidTimes[i],vidOnlyPeriod=True)}})
    return test
