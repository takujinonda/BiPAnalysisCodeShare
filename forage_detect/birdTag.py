from pathlib import Path
# bring in functions from main_func
# from forage_detect.main_func import *
# bring in functions from utils
import utils.analyseAcc as accFn
import utils.analyseGPS as gpsFn
import utils.DVLutils as dvlFn
import utils.loadIn as load
from math import pi
from statistics import median, mean
from typing import Union

import numpy as np
import pandas as pd

class birdTag:
    """
    Class for avian biologging tags. Current support includes AxyTrek and DVL with a look toward incorporating NinjaScan data. Requires file path(s) and string descriptor of tag ('Axy','DVL'). Methods include reading of data from raw CSV/txt or BiP database. Data loads can include GPS or not (where present). Vectorised haversine distance and speed calculation from GPS. Dynamic/static acceleration calculation via equiripple lowpass filter. Pitch calculation. Data removal from proximity to location.
    """
    def __init__(
            self,
            filepath: str,
            type: str,
            tagname: str,
            accfs: int,
            long_acc_name: Union[str,None] = None,
            gps_fixes_per_minute: Union[int,None] = None,
            accStart: Union[np.datetime64,None] = None,
            vidStart: Union[np.datetime64,None] = None
            ):
        self.filepath = filepath
        self.type = type.casefold()
        self.tagname = tagname
        self.accfs = accfs
        self.gps_fixes_per_minute = gps_fixes_per_minute
        self.long_acc_name = long_acc_name
        self.accStart = pd.to_datetime(accStart,format='%d/%m/%Y %H:%M:%S')
        if vidStart is not None:
            vidStart = pd.to_datetime(vidStart,format='%d/%m/%Y %H:%M:%S')
        self.vidStart = vidStart

    def readin(
            self,
            vidOnlyPeriod: bool = True
            ) -> None:
        if self.type == 'axy':
            self.gps = load.readAxy(self.filepath,cols='gps').reset_index()
            self.acc = load.readAxy(self.filepath,cols='acc').reset_index()
        if self.type == 'bip':
            self.gps = load.readBiP(self.filepath,col='gps').reset_index()
            self.acc = load.readBiP(self.filepath,col='acc').reset_index()
        if self.type == 'dvl':
            self.acc = load.readDVL(self.filepath, 
                                    accStart=self.accStart,
                                    fs=self.accfs,
                                    vidStart=self.vidStart,vidOnlyPeriod=vidOnlyPeriod).reset_index()

    def readBeh(
            self,
            behavPath: str
            ) -> None:
        self.dvl_beh = dvlFn.readBeh(behavPath,self.tagname)

    def accFeatures(
            self,
            passband: float = 1.5,
            stopband: float = 3.0
            ) -> None:
        """
        """
        # check if data read in
        if not hasattr(self, 'acc'):
            print("No acceleration data present")
        else:
            out = accFn.accFeatures(self.acc,self.long_acc_name,passband,stopband,self.accfs)
            out.index = self.acc.index
            self.acc = pd.concat([self.acc,out],axis=1)

    def rollSum(
            self,
            minfreq: float = 3.0,
            maxfreq: float = 5.0
            ) -> None:
        """
        """
        self.rolling_freq_sum = accFn.rollingSpecSum(self.acc.Z,minFreq=minfreq,maxFreq=maxfreq,fs=self.accfs,dur=60,inclusive=False)

    def flight_est(
            self,
            numPoints: int,
            removeErr: bool = False
            ) -> None: 
        # check if rollSum calculation has been performed
        if not hasattr(self, 'rolling_freq_sum'):
            print("Performing spectral rolling sum")
            self.rollSum()
        
        # if dvl tag, reduce magnitude of erroneous category
        if (self.type == 'dvl') & removeErr:
            self.flInds, self.flight = accFn.flightestimate(signal=self.acc.Z,rollSum=self.rolling_freq_sum,fs=self.accfs,behav_data=self.dvl_beh,dt=self.acc.DT,removeErr=removeErr,numPoints=numPoints)
        else:
            self.flInds, self.flight = accFn.flightestimate(signal=self.acc.Z,rollSum=self.rolling_freq_sum,fs=self.accfs,dt=self.acc.DT,numPoints=numPoints)

    def flight_thresholds(
            self
            ) -> None:
        if (not hasattr(self, 'flInds')) | (not hasattr(self, 'acc')):
            print(f"Flight estimate or acceleration data required")
        else:
            self.acc = pd.concat([self.acc,accFn.flight_est_thresholds(self.acc,self.flInds)],axis=1)

    def flapping(self):
        """
        Calculate mask of flapping behaviour and 'bouts', grouped by `flap_freq` and `bout_gap` seconds, respectively. If the tag is DVL, flapping peaks/troughs will be found within estimate flight periods. This function requires a flight mask to be present.
        """

        if not hasattr(self, 'flight') & (self.type == 'dvl'):
            print("Add estimated flight periods")
        
        if self.type == 'dvl':
            self.flap, self.flap_bouts, self.flap_start, self.flap_end = accFn.flap(sig = self.acc.Z, fs = self.accfs, bout_gap = 10, flap_freq = 4,find_in_flight_periods=True,flinds=self.flInds)
        else:
            self.flap, self.flap_bouts, self.flap_start, self.flap_end = accFn.flap(sig = self.acc.Z, fs = self.accfs, bout_gap = 10, flap_freq = 4)

    def pitchPT(
            self
            ) -> None:
        """
        Calculate pitch thresholds from estimated flight periods
        """
        if not hasattr(self, 'flInds'):
            print("Add estimated flight periods")
        
        self.pitFL = 1.5*accFn.flight_pitch_changes(self.acc.pitch,self.flInds,findVal='max')
        self.ODmFL = accFn.flight_pitch_changes(self.acc.ODmn,self.flInds,findVal='min')

    def beh_detect(
            self,
            toEx : Union[float,None]
            ) -> None:
        """
        """
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

        if (toEx is None) & (self.type == 'axy'):
            raise ValueError('Pitch median threshold not given')
        
        # define an empty ethogram
        self.EthBeh = ["" for _ in range(len(self.acc))]
        ODlow = np.where(self.ODmn < .2)[0]
        self.EthBeh[ODlow] = "Rest"
        DiveUp = median([np.mean(self.acc.pitch[x]) for x in self.flInds]) + 2*median([np.var(self.acc.pitch[x]) for x in self.flInds])
        DiveDown = median([np.min(self.acc.pitch[x]) for x in self.flInds]) - 30
        self.EthBeh[np.where(self.flap_bouts == 1)[0]] = 'FL'
        # find pitch changes
        pitpeaks,pittroughs,_ = accFn.peak_trough(self.acc.pitch)
        PitUp = self.pitch[pitpeaks] - self.pitch[pittroughs]
        PitDown = self.pitch[pittroughs] - self.pitch[pitpeaks]
        PitLarge = list(set(pittroughs[PitUp > toEx] + pitpeaks[PitDown > toEx]))
        PitUpL = pittroughs[PitUp > toEx]
        PitDownL = pitpeaks[PitDown > toEx]

        median_mean_flight_pitch = median([mean(self.pitch[x]) for x in self.flInds])
        # examine flight period starts and search for large pitch change
        for fl_start,fl_end in zip(self.flap_start,self.flap_end):
            if any((PitUpL > (max([fl_start - self.accfs*4,1]))) & (PitUpL < min([max(fl_start) + self.accfs*4,len(self.acc)]))):
                # find where the nearest increase in pitch is within 4 seconds of flight starting
                TkoStart = min([PitUpL[((PitUpL > (max([fl_start - self.accfs*4,1]))) & (PitUpL < min([max(fl_start) + self.accfs*4,len(self.acc)])))]])
                # find pitch peak after this increase in pitch
                if any(self.pitch[min(pitpeaks[pitpeaks > TkoStart]):] <= median_mean_flight_pitch):
                    TkoEnd = min([fl_end[fl_end > min(pitpeaks[pitpeaks > TkoStart]) + np.where(self.pitch[min(pitpeaks[pitpeaks > TkoStart]):] < median_mean_flight_pitch)[0][0]],TkoStart + self.accfs*5])
                    self.EthBeh[TkoStart:TkoEnd] = 'Takeoff'
        
        # identify foraging
        Pitdif = np.where(np.diff(PitLarge) > (23.3 * self.accfs))[0]
        PitOutSd = np.zeros(len(self.acc.pitch),dtype=int).tolist()
        # combine all foraging within 23.3 seconds
        for b in range(len(Pitdif) - 1):
            if b == 0:
                if PitLarge[Pitdif[b]] - PitLarge[0] > self.accfs*1.5:
                    PitOutSd[PitLarge[0]:PitLarge[Pitdif[0]]] = [1] * (PitLarge[Pitdif[0]] - PitLarge[0])
            elif PitLarge[Pitdif[b + 1]] - PitLarge[Pitdif[b] + 1] > self.accfs * 1.5:
                


# for DVL tags, need to calculate flight periods and thresholds across all tags to save to object. These can then be used for threshold analysis of individual tags. Define method for working across all tags, including reading data from all
# as we will be working with all for defining accuracies, this can be it's own dedicated part of the class. Can later separate etc, but for now just keep together

    def calculate_thresholds(
            self
            ) -> None:
        """
        Calculate all the required threshold for behaviour estimation.
        """
        # check if acceleration characteristics are present
        if not 'pitch' in self.acc:
            self.accFeatures()
        # check if flight estimate periods are present
        if not hasattr(self, 'flInds'):
            print('Estimating flight periods')
            if self.type == 'dvl':
                self.flight_est(numPoints=10,removeErr=True)
        # check if flight estimate median pitch range recorded
        if not hasattr(self, 'pitFL'):
            self.pitchPT()
        return 1.5 * self.pitFL