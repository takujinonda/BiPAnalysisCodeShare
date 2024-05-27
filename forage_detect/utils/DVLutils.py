import pandas as pd
import numpy as np

def readBeh(behaviour_file):
    out = pd.read_csv(behaviour_file, sep = ',', usecols = [0,1,2,4], dtype = {'Tag': str})
    # replace 'Dive' behaviour with 's' (surface) or 'd' (dive)
    out.loc[out.Behaviour == 'Dive','Behaviour'] = out.loc[out.Behaviour == 'Dive','ForageBeh']
    # remove superfluous column
    out.drop('ForageBeh',axis=1,inplace=True)

    # set time to datetime
    out.Time = pd.to_datetime('31/08/2018 ' + out.Time, format = '%d/%m/%Y %H:%M:%S.%f')

    return out