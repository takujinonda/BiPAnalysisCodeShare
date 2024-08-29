# Script to assess accuracy of behaviour detection method on DVL tag data

import platform 
import pandas as pd

import os
os.chdir('forage_detect')
from birdTag import birdTag

import glob

# parent folder for all DVL data
if platform.system() == 'Darwin':
    dvlFolders = "/Users/aran/Documents/DVL/"
else:
    dvlFolders = "H:/My Drive/Data/DVL/"

# read in the video behaviour records
if platform.system() == 'Darwin':
    behav_class_loc = '/Users/aran/Library/CloudStorage/GoogleDrive-aran.garrod@googlemail.com/My Drive/Data/DVL/DiveBehavClass.csv'
else:
    behav_class_loc = 'H:/My Drive/Data/DVL/DiveBehavClass.csv'

accFiles = glob.glob(dvlFolders + '*/*acc*.txt')

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
    # calculate flight estimates
    dvls[ind]['data'].flight_est(10,True)
    # calculate flight pitch changes
    dvls[ind]['data'].calculate_thresholds()
    dvls[ind]['data'].flapping()

# calculate pitch differences for all tags
medianPitchChanges = []
for x in dvls:
    medianPitchChanges.append(x['data'].calculate_thresholds())
from statistics import median
toEx = median(medianPitchChanges)

dvls[0]['data'].beh_detect(toEx)

import utils.analyseAcc as accFn

import numpy as np
from itertools import compress


dvls[0]['data'].EthBeh = ["" for _ in range(len(dvls[0]['data'].acc))]
ODlow = np.where(dvls[0]['data'].acc.ODmn < .2)[0]
for b in ODlow:
    dvls[0]['data'].EthBeh[b] = "Rest"
DiveUp = median([np.mean(dvls[0]['data'].acc.pitch[x]) for x in dvls[0]['data'].flInds]) + 2*median([np.var(dvls[0]['data'].acc.pitch[x]) for x in dvls[0]['data'].flInds])
DiveDown = median([np.min(dvls[0]['data'].acc.pitch[x]) for x in dvls[0]['data'].flInds]) - 30
for b in np.where(dvls[0]['data'].flap_bouts == 1)[0]:
    dvls[0]['data'].EthBeh[b] = "FL"
# find pitch changes
pitpeaks,pittroughs,_ = accFn.peak_trough(dvls[0]['data'].acc.pitch)
PitUp = np.array(dvls[0]['data'].acc.pitch[pitpeaks]) - np.array(dvls[0]['data'].acc.pitch[pittroughs])
PitDown = np.array(dvls[0]['data'].acc.pitch[pittroughs]) - np.array(dvls[0]['data'].acc.pitch[pitpeaks])
PitUpL = list(compress(pitpeaks,abs(PitDown) > toEx))
PitDownL = list(compress(pittroughs,PitUp > toEx))
PitLarge = sorted(PitUpL + PitUpL)

PitLarge = list(set(pittroughs[PitUp > toEx] + pitpeaks[PitDown > toEx]))
PitUpL = pittroughs[PitUp > toEx]
PitDownL = pitpeaks[PitDown > toEx]


[x for x in PitUp if str(x) != 'nan']

test = dvls[0]['data']

np.sign(np.array(test.flap_bouts)).diff(1).eq(1) & np.sign(np.array(test.flap_bouts).eq(1))
df[np.sign(df[0]).diff(1).eq(1) & np.sign(df[0]).eq(1)]


transition_points = []
transition_points.append(np.where(test.flap_bouts == 1)[0][0])
# transition back to 0
next = np.where(test.flap_bouts[transition_points[-1]:] == 0)[0][0]
test.flap_bouts[transition_points[-1]:transition_points[-1] + 28]


min([fl_end[fl_end > min(dvls[0].pitpeaks[dvls[0].pitpeaks > TKoStart]) + np.where(dvls[0].pitch[min(dvls[0].pitpeaks[dvls[0].pitpeaks > TKoStart]):] < median_mean_flight_pitch)[0][0]],TKoStart + dvls[0].accfs*5])


test.flap_end[test.flap_end > min(test.pitpeaks)]

np.where(test.flap_bouts[(transition_points[-1]+next):] == 1)[0][0]

def findTransitionPoint(arr): 
    transition_points = []
    # define the first point
    transition_points.append(np.where(arr == 1)[0][0])
    continue_searching = True
    while sum(np.diff(arr[transition_points[-1]:]) != 0):
        next = np.where(arr[transition_points[-1]:] == 0)[0][0]
        transition_points.append(np.where(arr[(transition_points[-1]+next):] == 1)[0][0])
    return transition_points

findTransitionPoint(test.flap_bouts)


    # Initialise lower and upper 
    # bounds 
    lb = 0
    ub = n - 1
  
    # Perform Binary search 
    while (lb <= ub): 
        # Find mid 
        mid = (int)((lb + ub) / 2) 
  
        # update lower_bound if 
        # mid contains 0 
        if (arr[mid] == 0): 
            lb = mid + 1
  
        # If mid contains 1 
        elif (arr[mid] == 1): 
              
            # Check if it is the  
            # left most 1 Return 
            # mid, if yes 
            if (mid == 0 or (mid > 0 and arr[mid - 1] == 0)): 
                return mid 
  
            # Else update  
            # upper_bound 
            ub = mid-1
      
    return -1
  
# Driver code 
arr = [0, 0, 0, 0, 1, 1] 
n = len(arr) 
point = findTransitionPoint(test.flap_bouts, len(test.flap_bouts)); 
if(point >= 0): 
    print("Transition point is ", point) 
else: 
    print("There is no transition point") 
  