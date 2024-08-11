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
median(medianPitchChanges)


test = dvls[0]['data']

np.sign(np.array(test.flap_bouts)).diff(1).eq(1) & np.sign(np.array(test.flap_bouts).eq(1))
df[np.sign(df[0]).diff(1).eq(1) & np.sign(df[0]).eq(1)]


transition_points = []
transition_points.append(np.where(test.flap_bouts == 1)[0][0])
# transition back to 0
next = np.where(test.flap_bouts[transition_points[-1]:] == 0)[0][0]
test.flap_bouts[transition_points[-1]:transition_points[-1] + 28]


min([fl_end[fl_end > min(self.pitpeaks[self.pitpeaks > TKoStart]) + np.where(self.pitch[min(self.pitpeaks[self.pitpeaks > TKoStart]):] < median_mean_flight_pitch)[0][0]],TKoStart + self.accfs*5])


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
  