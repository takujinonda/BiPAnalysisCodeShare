from pathlib import Path
# bring in functions from main_func
# from forage_detect.main_func import *
# bring in functions from utils
import utils.analyseAcc as accFn
import utils.analyseGPS as gpsFn
import utils.DVLutils as dvlFn
import utils.loadIn as load
from math import pi
from statistics import median

import numpy as np
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
    def __init__(self,filepath,type,tagname,accfs,long_acc_name=None,gpsfpm=None,accStart=None,vidStart=None):
        self.filepath = filepath
        self.type = type
        self.tagname = tagname
        self.accfs = accfs
        self.gpsfpm = gpsfpm
        self.long_acc_name = long_acc_name
        self.accStart = pd.to_datetime(accStart,format='%d/%m/%Y %H:%M:%S')
        self.vidStart = pd.to_datetime(vidStart,format='%d/%m/%Y %H:%M:%S')

    def readin(self,vidOnlyPeriod=True):
        if self.type.lower() == 'axy':
            self.gps = load.readAxy(self.filepath,cols='gps').reset_index()
            self.acc = load.readAxy(self.filepath,cols='acc').reset_index()
        elif self.type.lower() == 'bip':
            self.gps = load.readBiP(self.filepath,col='gps').reset_index()
            self.acc = load.readBiP(self.filepath,col='acc').reset_index()
        elif self.type.lower() == 'dvl':
            self.acc = load.readDVL(self.filepath, accStart=self.accStart, fs=self.accfs, vidStart=self.vidStart,vidOnlyPeriod=vidOnlyPeriod).reset_index()

    def readBeh(self,behavPath):
        self.dvl_beh = dvlFn.readBeh(behavPath,self.tagname)

    def accFeatures(self,passband=1.5,stopband=3):
        # check if data read in
        if not hasattr(self, 'acc'):
            print("No acceleration data present")
        else:
            out = accFn.accFeatures(self.acc,self.long_acc_name,passband,stopband,self.accfs).reset_index()
            self.acc = pd.concat([self.acc,out],axis=1)

    def rollSum(self,minfreq=3,maxfreq=5):
        self.rolling_freq_sum = accFn.rollingSpecSum(self.acc.Z,minFreq=minfreq,maxFreq=maxfreq,fs=self.accfs,dur=60,inclusive=False)

    def flight_est(self,numPoints,removeErr=False):
        # check if rollSum calculation has been performed
        if not hasattr(self, 'rolling_freq_sum'):
            print("Performing spectral rolling sum")
            self.rolling_freq_sum()
        
        # if dvl tag, reduce magnitude of erroneous category
        if (self.type.lower() == 'dvl') & removeErr:
            self.flInds, self.flight = accFn.flightestimate(signal=self.acc.Z,rollSum=self.rolling_freq_sum,fs=self.accfs,behav_data=self.dvl_beh,dt=self.acc.DT,removeErr=removeErr,numPoints=numPoints)
        else:
            self.flInds, self.flight = accFn.flightestimate(signal=self.acc.Z,rollSum=self.rolling_freq_sum,fs=self.accfs,dt=self.acc.DT,numPoints=numPoints)

    def flight_thresholds(self):
        if (not hasattr(self, 'flInds')) | (not hasattr(self, 'acc')):
            print(f"Flight estimate or acceleration data required")
        else:
            self.acc = pd.concat([self.acc,accFn.flight_est_thresholds(self.acc,self.flInds)],axis=1)

    def flapping(self):
        """Calculate mask of flapping behaviour and 'bouts', grouped by `flap_freq` and `bout_gap` seconds, respectively. If the tag is DVL, flapping peaks/troughs will be found within estimate flight periods. This function requires a flight mask to be present.
        """

        if not hasattr(self, 'flight') & (self.type.lower() == 'dvl'):
            print("Add estimated flight periods")
        
        if self.type.lower() == 'dvl':
            self.flap, self.flap_bouts = accFn.flap(sig = self.acc.Z, fs = self.accfs, bout_gap = 10, flap_freq = 4,find_in_flight_periods=True,flinds=self.flInds)
        else:
            self.flap, self.flap_bouts = accFn.flap(sig = self.acc.Z, fs = self.accfs, bout_gap = 10, flap_freq = 4)

    def pitchPT(self):
        """Calculate pitch thresholds from estimated flight periods
        """
        if not hasattr(self, 'flInds'):
            print("Add estimated flight periods")
        
        self.pitFL = 1.5*accFn.flight_pitch_changes(self.acc.pitch,self.flInds,findVal='max')
        self.ODmFL = accFn.flight_pitch_changes(self.acc.ODmn,self.flInds,findVal='min')

    def beh_detect(self):#
        if not hasattr(self, 'pitFL'):
            print('Pre-requesite calculations not performed\nRunning now...')
            if not hasattr(self, 'acc'):
                self.readin()
                self.accFeatures()
            if not hasattr(self, 'rolling_freq_sum'):
                print('Rolling spectral sum not calculated\nRunning now...')
                self.rollSum()
            if not hasattr(self, 'flight'):
                print('Flight not calculated\nRunning now...')
                self.flight_est()
            if not hasattr(self,'ODmFL') & (self.type == 'dvl'):
                print('Calculating some thresholds')
                self.pitchPT()
        
        # define an empty ethogram
        self.EthBeh = ["" for x in range(len(self.acc))]
        ODlow = np.where(self.ODmn < .2)[0]
        self.EthBeh[ODlow] = "Rest"
        DiveUp = median([np.mean(self.acc.pitch[x]) for x in self.flInds]) + 2*median([np.var(self.acc.pitch[x]) for x in self.flInds])
        DiveDown = median([np.min(self.acc.pitch[x]) for x in self.flInds]) - 30
        self.EthBeh[np.where(self.flap_bouts == 1)[0]] = 'FL'
        # find pitch changes
        pitpeaks,pittroughs,_ = accFn.peak_trough(self.acc.pitch)
        # change to start with trough
        pitpeaks = pitpeaks[1:]
        pittroughs = pittroughs[:-2]
        data = self.acc.pitch[pitpeaks].values - self.acc.pitch[pittroughs].values
        PitLarge = data > self.pitFL
        pitpeaks = pitpeaks[data > self.pitFL]

# for DVL tags, need to calculate flight periods and thresholds across all tags to save to object. These can then be used for threshold analysis of individual tags. Define method for working across all tags, including reading data from all
# as we will be working with all for defining accuracies, this can be it's own dedicated part of the class. Can later separate etc, but for now just keep together