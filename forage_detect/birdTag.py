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
from operator import itemgetter
from itertools import compress

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

    @staticmethod
    def find_changes(
        list_to_search: list,
        value_to_change
        ) -> list[list,list]:
        """
        Find the start and end points of consecutive series within
        `list_to_search` list. Starts and ends should be indeces of where
        `list_to_search` is equal to `value_to_change`.
        """
        # check if value to change is present in list
        if not set([value_to_change]).issubset(set(list_to_search)):
            return None, None
        is_desired = [x == value_to_change for x in list_to_search]
        starts = []
        ends = []
        if list_to_search[0] == value_to_change:
            starts.insert(0,0)
        [starts.append(i+1) for i, x in enumerate([(sub2 - sub1) == 1 for sub1,sub2 in zip(is_desired, is_desired[1:])]) if x]
        [ends.append(i+1) for i, x in enumerate([(sub1 - sub2) == 1 for sub1,sub2 in zip(is_desired, is_desired[1:])]) if x]
        if list_to_search[-1] == value_to_change:
            ends.append(len(list_to_search))

        return starts, ends

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
        ODlow = np.where(self.acc.ODmn < .2)[0]
        for b in ODlow:
            self.EthBeh[b] = "Rest"
        DiveUp = median([np.mean(self.acc.pitch[x]) for x in self.flInds]) + 2*median([np.var(self.acc.pitch[x]) for x in self.flInds])
        DiveDown = median([np.min(self.acc.pitch[x]) for x in self.flInds]) - 30
        for b in np.where(self.flap_bouts == 1)[0]:
            self.EthBeh[b] = "FL"
        # find pitch changes
        pitpeaks,pittroughs,_ = accFn.peak_trough(self.acc.pitch)
        PitUp = np.array(self.acc.pitch[pitpeaks]) - np.array(self.acc.pitch[pittroughs])
        PitDown = np.array(self.acc.pitch[pittroughs]) - np.array(self.acc.pitch[pitpeaks])
        PitUpL = list(compress(pitpeaks,abs(PitDown) > toEx))
        PitDownL = list(compress(pittroughs,PitUp > toEx))
        PitLarge = sorted(PitUpL + PitUpL)

        median_mean_flight_pitch = median([np.mean(self.acc.pitch[x]) for x in self.flInds])
        # examine flight period starts and search for large pitch change
        for fl_start,fl_end in zip(self.flap_start,self.flap_end):
            idx_range = range(max([fl_start - self.accfs*4,1]),min([fl_start + self.accfs*4,len(self.acc)]))
            if any(list(filter(lambda a: a in set(PitUpL),set(list(idx_range))))):
                # find where the nearest increase in pitch is within 4 seconds of flight starting
                TkoStart = min([PitUpL[((PitUpL > (max([fl_start - self.accfs*4,1]))) & (PitUpL < min([max(fl_start) + self.accfs*4,len(self.acc)])))]])
                # find pitch peak after this increase in pitch
                if any(self.acc.pitch[min(pitpeaks[pitpeaks > TkoStart]):] <= median_mean_flight_pitch):
                    TkoEnd = min([fl_end[fl_end > min(pitpeaks[pitpeaks > TkoStart]) + np.where(self.acc.pitch[min(pitpeaks[pitpeaks > TkoStart]):] < median_mean_flight_pitch)[0][0]],TkoStart + self.accfs*5])
                    self.EthBeh[TkoStart:TkoEnd] = ['Takeoff'] * (TkoEnd-TkoStart)
        
        # identify foraging
        Pitdif = np.where(np.diff(PitLarge) > (23.3 * self.accfs))[0]
        PitOutSd = np.zeros(len(self.acc.pitch),dtype=int).tolist()
        # combine all foraging within 23.3 seconds
        for b in range(len(Pitdif) - 1):
            if b == 0:
                if PitLarge[Pitdif[b]] - PitLarge[0] > self.accfs*1.5:
                    PitOutSd[PitLarge[0]:PitLarge[Pitdif[0]]] = [1] * (PitLarge[Pitdif[0]] - PitLarge[0])
            elif PitLarge[Pitdif[b + 1]] - PitLarge[Pitdif[b] + 1] > self.accfs * 1.5:
                PitOutSd[PitLarge[Pitdif[b] + 1]:PitLarge[Pitdif[b+1]]] = [1] * (PitLarge[Pitdif[b+1]] - PitLarge[Pitdif[b] + 1])
                if b == (len(Pitdif) - 2):
                    PitOutSd[PitLarge[Pitdif[-1]+1]:PitLarge[-1]] = [1] * (PitLarge[-1] - PitLarge[Pitdif[-1]])

        all_pass = [sum(x) for x in zip(PitOutSd, [x == 'FL' for x in self.EthBeh], [x > self.ODmFL for x in self.acc.ODmn])]
        all_pass_indeces = [i for i, x in enumerate(all_pass) if x == 3]
        for x in all_pass_indeces:
            self.EthBeh[x] = "Forage"

        dives = np.zeros(len(self.acc.pitch),dtype=int).tolist()
        ForSt, ForEd = self.find_changes(self.EthBeh, 'Forage')
        # dives will have a significant drop in pitch to start and are followed
        # by an increase later as the bird returns to the surface
        for fs,fe in zip(ForSt,ForEd):
            # search for dives prior to foraging indication
            if any(self.acc.pitch[max(1,(fs-self.accfs)):fe] < DiveDown):
                dive_from = next((i for i in reversed(pittroughs) if i < fs))
                # check the dive has a large downward pitch to remove possible
                # take-offs
                if not set([next((i for i in reversed(pitpeaks) if i < dive_from))]).issubset(set(PitDownL)):
                    continue
                else:
                    # search if pitch following potential dive reaches above
                    # normal levels
                    if any(self.acc.pitch[dive_from:(dive_from + (self.accfs*10))] > DiveUp):
                        dive_to = pitpeaks[np.argmax((pitpeaks > dive_from) & (self.acc.pitch[pitpeaks] > DiveUp))]
                        dives[dive_from:dive_to] = [1]*(dive_to - dive_from)
            if any(self.acc.pitch[fs:fe] < DiveDown):
                troughsin = pittroughs[((pittroughs >= fs) & (pittroughs <= fe))]
                diveTroughs = troughsin[self.acc.pitch[troughsin] < DiveDown]
                for dt in diveTroughs:
                    if any(self.acc.pitch[dt:(dt + (self.accfs * 10))] > DiveUp):
                        if not set([next((i for i in reversed(pitpeaks) if i < dt))]).issubset(set(PitDownL)):
                            continue
                        else:
                            dive_to = pitpeaks[np.argmax((pitpeaks > dt) & (self.acc.pitch[pitpeaks] > DiveUp))]
                            dives[dt:dive_to] = [1] * (dive_to - dt)
        for d in np.where(dives == 1)[0]:
            self.EthBeh[d] = "Dive"
            
        # remove foraging bouts where fewer than two large downward pitch
        # changes occur within 1s of each other
        if ForSt != []:
            for fs,fe in zip(ForSt, ForEd):
                if sum(np.diff(PitUpL[(PitUpL >= fs) & (PitUpL <= fe)]) < (self.accfs * 2)) < 2:
                    self.EthBeh[fs:fe] = ["Failed 1"] * (fe - fs)
                if sum(np.diff(PitDownL[(PitDownL >= fs) & (PitDownL <= fe)]) < (self.accfs * 2)) < 2:
                    self.EthBeh[fs:fe] = ["Failed 1"] * (fe - fs)
        

    # for DVL tags, need to calculate flight periods and thresholds across all
    # tags to save to object. These can then be used for threshold analysis of
    # individual tags. Define method for working across all tags, including
    # reading data from all as we will be working with all for defining
    # accuracies, this can be it's own dedicated part of the class. Can later
    # separate etc, but for now just keep together

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