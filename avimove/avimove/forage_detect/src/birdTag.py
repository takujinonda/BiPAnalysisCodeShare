import avimove.forage_detect.src.utils.analyseAcc as accFn
import avimove.forage_detect.src.utils.analyseGPS as gpsFn
import avimove.forage_detect.src.utils.DVLutils as dvlFn
import avimove.forage_detect.src.utils.loadIn as load
import matplotlib.pyplot as plt
import distinctipy
import numpy as np
import pandas as pd

# suppress warnings on chained assignment (not relevant here)
pd.options.mode.chained_assignment = None
import datetime as dt
from typing import Union

from math import pi
from statistics import median, mean
from typing import Union
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
from itertools import compress


class birdTag:
    """
    Class for avian biologging tags. Current support includes AxyTrek and DVL
    with a look toward incorporating NinjaScan data. Requires file path(s) and
    string descriptor of tag ('Axy','DVL'). Methods include reading of data from
    raw CSV/txt or BiP database. Data loads can include GPS or not (where
    present). Vectorised haversine distance and speed calculation from GPS.
    Dynamic/static acceleration calculation via equiripple lowpass filter. Pitch
    calculation. Data removal from proximity to location.
    """

    def __init__(
        self,
        filepath: str,
        tag_type: str,
        tagname: str,
        accfs: int,
        long_acc_name: Union[str, None] = None,
        gps_fixes_per_minute: Union[int, None] = None,
        accStart: Union[np.datetime64, None] = None,
        vidStart: Union[np.datetime64, None] = None,
        *args,
        **kwargs,
    ):
        """
        Args
        ----
        filepath
            Path to signal data
        tag_type
            Type of tag. Currently choice of 'DVL' and 'AxyTrek' (case
            insensitive)
        tagname
            Tag identifier
        accfs
            Acceleration sampling frequency (Hz).
        long_acc_name
            list of names for acceleration signals. Must conform to the
            following order - longitudinal (along body axis), dorsoventral
            (vertical axis), lateral
        gps_fixes_per_minute
            Number of expected GPS fixes per minute
        accStart
            DateTime start of acceleration signal recording
        vidStart
            DateTime start of video recording
        """
        self.filepath = filepath
        self.tag_type = tag_type.casefold()
        self.tagname = tagname
        self.accfs = accfs
        self.gps_fixes_per_minute = gps_fixes_per_minute
        self.long_acc_name = long_acc_name
        self.accStart = pd.to_datetime(accStart, format="%d/%m/%Y %H:%M:%S")
        if vidStart is not None:
            vidStart = pd.to_datetime(vidStart, format="%d/%m/%Y %H:%M:%S")
        self.vidStart = vidStart

    @staticmethod
    def find_changes(list_to_search: list, value_to_change) -> list[list, list]:
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
            starts.insert(0, 0)
        [
            starts.append(i + 1)
            for i, x in enumerate(
                [(sub2 - sub1) == 1 for sub1, sub2 in zip(is_desired, is_desired[1:])]
            )
            if x
        ]
        [
            ends.append(i + 1)
            for i, x in enumerate(
                [(sub2 - sub1) == -1 for sub1, sub2 in zip(is_desired, is_desired[1:])]
            )
            if x
        ]
        if list_to_search[-1] == value_to_change:
            ends.append(len(list_to_search))

        return starts, ends

    def readin(self, vidOnlyPeriod: bool = True) -> None:
        if self.tag_type == "axy":
            self.gps = load.readAxy(self.filepath, cols="gps").reset_index()
            self.acc = load.readAxy(self.filepath, cols="acc").reset_index()
        if self.tag_type == "bip":
            self.gps = load.readBiP(self.filepath, col="gps").reset_index()
            self.acc = load.readBiP(self.filepath, col="acc").reset_index()
        if self.tag_type == "dvl":
            self.acc = load.readDVL(
                self.filepath,
                accStart=self.accStart,
                fs=self.accfs,
                vidStart=self.vidStart,
                vidOnlyPeriod=vidOnlyPeriod,
            ).reset_index()

    def readBeh(self, behavPath: str) -> None:
        self.dvl_beh = dvlFn.readBeh(behavPath, self.tagname)

    @staticmethod
    def round_seconds(obj: dt.datetime) -> dt.datetime:
        """
        Round a datetime to the nearest second.
        """
        if obj.microsecond >= 500_000:
            obj += dt.timedelta(seconds=1)
        return obj.replace(microsecond=0)

    def create_regularised_behaviours(self):
        """
        Create a regularly sampled time series of video-recorded behaviours. Simply
        assign the last behaviour before the second ends.
        """
        if self.tag_type != "dvl":
            pass
        else:
            if not hasattr(self, "dvl_beh"):
                print("Behaviour data not read in")
            self.dvl_beh["rounded_time"] = [
                self.round_seconds(x) for x in self.dvl_beh.Time
            ]

            start_time = self.dvl_beh.Time.iloc[0]
            end_time = self.dvl_beh.Time.iloc[0] + dt.timedelta(hours=2)
            reg_date_range = pd.date_range(start=start_time, end=end_time, freq="1s")
            reg_beh = np.repeat("not set", len(reg_date_range))
            from_ind = 0
            for idx, x in enumerate(self.dvl_beh.rounded_time[:-2]):
                to_ind = np.argmax(reg_date_range > x) - 1
                reg_beh[from_ind:to_ind] = self.dvl_beh.Behaviour.iloc[idx]
                from_ind = to_ind
            reg_beh[from_ind:] = self.dvl_beh.Behaviour.iloc[-1]

            regularised_behaviour = pd.DataFrame(
                {"time": reg_date_range, "beh": reg_beh}
            )

            self.regularised_behaviour = regularised_behaviour

    def upsample_behaviours(self, behaviours_duration: int = 2) -> None:
        """
        Upsample behaviour records to sampling frequency of acceleration records.
        Args
        ----
        behaviours_duration
            Duration (in hours) of behaviour records. Defaults to 2 for DVL.
        """
        # check behaviour records exist
        if not hasattr(self, "dvl_beh"):
            raise AttributeError
            print("Tag behaviours have not been read in")

        upsampled = self.dvl_beh[["Behaviour", "Time"]]
        upsampled.loc[-1] = [
            upsampled.Behaviour.iloc[-1],
            upsampled.Time[0]
            + pd.Timedelta(2, "hour")
            - pd.Timedelta(1 / self.accfs, "s"),
        ]
        upsampled.reset_index(drop=True)
        upsampled = upsampled.resample(pd.Timedelta(50, "ms"), on="Time").last().ffill()

        self.upsampled_beh = upsampled

    def accFeatures(self, passband: float = 1.5, stopband: float = 3.0) -> None:
        """ """
        # check if data read in
        if not hasattr(self, "acc"):
            print("No acceleration data present")
        else:
            out = accFn.accFeatures(
                self.acc, self.long_acc_name, passband, stopband, self.accfs
            )
            out.index = self.acc.index
            self.acc = pd.concat([self.acc, out], axis=1)

    def rollSum(self, minfreq: float = 3.0, maxfreq: float = 5.0) -> None:
        """ """
        self.rolling_freq_sum = accFn.rollingSpecSum(
            self.acc.Z,
            minFreq=minfreq,
            maxFreq=maxfreq,
            fs=self.accfs,
            dur=60,
            inclusive=False,
        )

    def flight_est(self, numPoints: int, removeErr: bool = False) -> None:
        # check if rollSum calculation has been performed
        if not hasattr(self, "rolling_freq_sum"):
            print("Performing spectral rolling sum")
            self.rollSum()

        # if dvl tag, reduce magnitude of erroneous category
        if (self.tag_type == "dvl") & removeErr:
            self.flInds, self.flight = accFn.flightestimate(
                signal=self.acc.Z,
                rollSum=self.rolling_freq_sum,
                fs=self.accfs,
                behav_data=self.dvl_beh,
                dt=self.acc.DT,
                removeErr=removeErr,
                numPoints=numPoints,
            )
        else:
            self.flInds, self.flight = accFn.flightestimate(
                signal=self.acc.Z,
                rollSum=self.rolling_freq_sum,
                fs=self.accfs,
                dt=self.acc.DT,
                numPoints=numPoints,
            )

    def flight_thresholds(self) -> None:
        if (not hasattr(self, "flInds")) | (not hasattr(self, "acc")):
            print(f"Flight estimate or acceleration data required")
        else:
            self.acc = pd.concat(
                [self.acc, accFn.flight_est_thresholds(self.acc, self.flInds)], axis=1
            )

    @staticmethod
    def flatten(xss):
        return [x for xs in xss for x in xs]

    def flapping(self):
        """
        Calculate mask of flapping behaviour and 'bouts', grouped by `flap_freq` and `bout_gap` seconds, respectively. If the tag is DVL, flapping peaks/troughs will be found within estimate flight periods. This function requires a flight mask to be present.
        """

        if not hasattr(self, "flight") & (self.tag_type == "dvl"):
            print("Add estimated flight periods")

        if self.tag_type == "dvl":
            self.flap, self.flap_bouts, self.flap_start, self.flap_end = accFn.flap(
                sig=self.acc.Z,
                fs=self.accfs,
                bout_gap=10,
                flap_freq=4,
                find_in_flight_periods=True,
                flinds=self.flInds,
            )
        else:
            self.flap, self.flap_bouts, self.flap_start, self.flap_end = accFn.flap(
                sig=self.acc.Z, fs=self.accfs, bout_gap=10, flap_freq=4
            )

    def pitchPT(self) -> None:
        """
        Calculate pitch thresholds from estimated flight periods
        """
        if not hasattr(self, "flInds"):
            print("Add estimated flight periods")

        self.pitFL = 1.5 * accFn.flight_pitch_changes(
            self.acc.pitch, self.flInds, findVal="max"
        )
        self.ODmFL = accFn.flight_pitch_changes(
            self.acc.ODmn, self.flInds, findVal="min"
        )

    def beh_detect(self, toEx: Union[float, None]) -> None:
        """ """
        if not hasattr(self, "pitFL"):
            print("Pre-requesite calculations not performed\nRunning now...")
            if not hasattr(self, "acc"):
                self.readin()
                self.accFeatures()
            if not hasattr(self, "rolling_freq_sum"):
                print("Rolling spectral sum not calculated\nRunning now...")
                self.rollSum()
            if not hasattr(self, "flight"):
                print("Flight not calculated\nRunning now...")
                self.flight_est()
            if not hasattr(self, "ODmFL") & (self.tag_type == "dvl"):
                print("Calculating some thresholds")
                self.pitchPT()

        if (toEx is None) & (self.tag_type == "axy"):
            raise ValueError("Pitch median threshold not given")

        # define an empty ethogram
        self.EthBeh = np.array(["Unknown" for _ in range(len(self.acc))])
        self.EthBeh[np.where(self.acc.ODmn < 0.2)] = "Rest"
        DiveUp = median([np.mean(self.acc.pitch[x]) for x in self.flInds]) + 2 * median(
            [np.var(self.acc.pitch[x]) for x in self.flInds]
        )
        DiveDown = median([np.min(self.acc.pitch[x]) for x in self.flInds]) - 30
        self.EthBeh[np.where(self.flap_bouts == 1)] = "FL"
        # find pitch changes
        pitpeaks, pittroughs, _ = accFn.peak_trough(self.acc.pitch)
        PitUp = np.array(self.acc.pitch[pitpeaks]) - np.array(
            self.acc.pitch[pittroughs]
        )
        PitDown = np.array(self.acc.pitch[pittroughs]) - np.array(
            self.acc.pitch[pitpeaks]
        )
        PitUpL = list(compress(pitpeaks, abs(PitDown) > toEx))
        PitDownL = list(compress(pittroughs, PitUp > toEx))
        PitLarge = sorted(PitUpL + PitUpL)

        # DROP THE TAKE-OFF AS THESE WERE NOT DEFINED FROM VIDEO FOOTAGE
        # median_mean_flight_pitch = median([np.mean(self.acc.pitch[x]) for x in self.flInds])
        # # examine flight period starts and search for large pitch change
        # for fl_start,fl_end in zip(self.flap_start,self.flap_end):
        #     idx_range = range(max([fl_start - self.accfs*4,1]),min([fl_start + self.accfs*4,len(self.acc)]))
        #     if any(list(filter(lambda a: a in set(PitUpL),set(list(idx_range))))):
        #         # find where the nearest increase in pitch is within 4 seconds of flight starting
        #         TkoStart = min([PitUpL[((PitUpL > (max([fl_start - self.accfs*4,1]))) & (PitUpL < min([max(fl_start) + self.accfs*4,len(self.acc)])))]])
        #         # find pitch peak after this increase in pitch
        #         if any(self.acc.pitch[min(pitpeaks[pitpeaks > TkoStart]):] <= median_mean_flight_pitch):
        #             TkoEnd = min([fl_end[fl_end > min(pitpeaks[pitpeaks > TkoStart]) + np.where(self.acc.pitch[min(pitpeaks[pitpeaks > TkoStart]):] < median_mean_flight_pitch)[0][0]],TkoStart + self.accfs*5])
        #             self.EthBeh[TkoStart:TkoEnd] = ['Takeoff'] * (TkoEnd-TkoStart)

        # identify foraging
        Pitdif = np.where(np.diff(PitLarge) > (23.3 * self.accfs))[0]
        PitOutSd = np.zeros(len(self.acc.pitch), dtype=int).tolist()
        # combine all foraging within 23.3 seconds
        for b in range(len(Pitdif) - 1):
            if b == 0:
                if PitLarge[Pitdif[b]] - PitLarge[0] > self.accfs * 1.5:
                    PitOutSd[PitLarge[0] : PitLarge[Pitdif[0]]] = [1] * (
                        PitLarge[Pitdif[0]] - PitLarge[0]
                    )
            elif PitLarge[Pitdif[b + 1]] - PitLarge[Pitdif[b] + 1] > self.accfs * 1.5:
                PitOutSd[PitLarge[Pitdif[b] + 1] : PitLarge[Pitdif[b + 1]]] = [1] * (
                    PitLarge[Pitdif[b + 1]] - PitLarge[Pitdif[b] + 1]
                )
                if b == (len(Pitdif) - 2):
                    PitOutSd[PitLarge[Pitdif[-1] + 1] : PitLarge[-1]] = [1] * (
                        PitLarge[-1] - PitLarge[Pitdif[-1]]
                    )

        all_pass = [
            sum(x)
            for x in zip(
                PitOutSd,
                [x == "FL" for x in self.EthBeh],
                [x > self.ODmFL for x in self.acc.ODmn],
            )
        ]
        all_pass_indeces = [i for i, x in enumerate(all_pass) if x == 3]
        for x in all_pass_indeces:
            self.EthBeh[x] = "Forage"

        dives = np.zeros(len(self.acc.pitch), dtype=int).tolist()
        ForSt, ForEd = self.find_changes(self.EthBeh, "Forage")
        # dives will have a significant drop in pitch to start and are followed
        # by an increase later as the bird returns to the surface
        for fs, fe in zip(ForSt, ForEd):
            # search for dives prior to foraging indication
            if any(self.acc.pitch[max(1, (fs - self.accfs)) : fe] < DiveDown):
                dive_from = next((i for i in reversed(pittroughs) if i < fs))
                # check the dive has a large downward pitch to remove possible
                # take-offs
                if not set(
                    [next((i for i in reversed(pitpeaks) if i < dive_from))]
                ).issubset(set(PitDownL)):
                    continue
                else:
                    # search if pitch following potential dive reaches above
                    # normal levels
                    if any(
                        self.acc.pitch[dive_from : (dive_from + (self.accfs * 10))]
                        > DiveUp
                    ):
                        dive_to = pitpeaks[
                            np.argmax(
                                (pitpeaks > dive_from)
                                & (self.acc.pitch[pitpeaks] > DiveUp)
                            )
                        ]
                        dives[dive_from:dive_to] = [1] * (dive_to - dive_from)
            if any(self.acc.pitch[fs:fe] < DiveDown):
                troughsin = pittroughs[((pittroughs >= fs) & (pittroughs <= fe))]
                diveTroughs = troughsin[self.acc.pitch[troughsin] < DiveDown]
                for dt in diveTroughs:
                    if any(self.acc.pitch[dt : (dt + (self.accfs * 10))] > DiveUp):
                        if not set(
                            [next((i for i in reversed(pitpeaks) if i < dt))]
                        ).issubset(set(PitDownL)):
                            continue
                        else:
                            dive_to = pitpeaks[
                                np.argmax(
                                    (pitpeaks > dt)
                                    & (self.acc.pitch[pitpeaks] > DiveUp)
                                )
                            ]
                            dives[dt:dive_to] = [1] * (dive_to - dt)
        for d in np.where(dives == 1)[0]:
            self.EthBeh[d] = "Dive"

        # remove foraging bouts where fewer than two large downward pitch
        # changes occur within 1s of each other
        # if ForSt != []:
        #     for fs,fe in zip(ForSt, ForEd):
        #         if sum(np.diff([PitUpL[x] for x in np.where([(x >= fs) & (x <= fe) for x in PitUpL])[0]]) > (self.accfs * 2)) < 2:
        #             self.EthBeh[fs:fe] = ["Failed 1"] * (fe - fs)
        #         if sum(np.diff([PitDownL[x] for x in np.where([(x >= fs) & (x <= fe) for x in PitDownL])[0]]) > (self.accfs * 2)) < 2:
        #             self.EthBeh[fs:fe] = ["Failed 1"] * (fe - fs)

    # for DVL tags, need to calculate flight periods and thresholds across all
    # tags to save to object. These can then be used for threshold analysis of
    # individual tags. Define method for working across all tags, including
    # reading data from all as we will be working with all for defining
    # accuracies, this can be it's own dedicated part of the class. Can later
    # separate etc, but for now just keep together

    def calculate_thresholds(self, threshold_scale: float = 1.5) -> None:
        """
        Calculate all the required threshold for behaviour estimation.
        """
        # check if acceleration characteristics are present
        if not "pitch" in self.acc:
            self.accFeatures()
        # check if flight estimate periods are present
        if not hasattr(self, "flInds"):
            print("Estimating flight periods")
            if self.tag_type == "dvl":
                self.flight_est(numPoints=10, removeErr=True)
        # check if flight estimate median pitch range recorded
        if not hasattr(self, "pitFL"):
            self.pitchPT()
        return threshold_scale * self.pitFL

    @staticmethod
    def get_changes_in_string_list(
        string_list: list, line_list: Union[list, None] = None
    ):
        """
        Convert list of continuous sections (string or other) into start/end
        indeces and corresponding value.
        """
        if line_list is not None:
            if len(string_list) != len(line_list):
                raise ValueError("string_list and line_list are not of equal length")
            to_use = line_list
        else:
            to_use = string_list
        # ensure the string_list which you are separating by is ndarray
        if type(string_list) != np.ndarray:
            string_list = np.array(string_list)
        change_idx = np.where(string_list[:-1] != string_list[1:])
        change_idx = np.append(change_idx, len(string_list) - 1)
        change_strings = string_list[change_idx]

        return change_idx, change_strings

    def get_lines_from_string_list(
        self, string_list: list, line_list: Union[list, None] = None
    ):
        """
        Convert list of continuous sections (string or other) and convert to
        structure for LineCollection plotting.
        """
        change_idx, change_strings = self.get_changes_in_string_list(
            string_list=string_list, line_list=line_list
        )

        if line_list is not None:
            if len(string_list) != len(line_list):
                raise ValueError("string_list and line_list are not of equal length")
            to_use = line_list
        else:
            to_use = string_list
        strings = np.empty(len(to_use), dtype=np.array(string_list).dtype)
        lines = []
        strings = []
        inds = []
        start = 0
        for idx, string in zip(change_idx, change_strings):
            strings[start : idx + 1] = [string]
            inds.append(list(range(start, idx + 1)))
            lines.append(
                np.array(
                    [
                        (x, y)
                        for x, y in zip(
                            list(range(start, idx + 1)),
                            [line_list[b] for b in list(range(start, idx + 1))],
                        )
                    ]
                )
            )
            start = idx + 1

        return inds, lines, strings

    @staticmethod
    def make_proxy(color: list, **kwargs):
        """
        Make custom lines for plot legend.

        Args
        ----
        color
            List of rgb array(s).
        """
        return Line2D([0, 1], [0, 1], color=color, **kwargs)

    def plot_acc_behaviours(self, acc_sig, cols=None, plot_vid_forage: bool = True):

        # check for simplified upsampled behaviours and generate if not present
        if (not hasattr(self,'upsampled_beh')) and (self.tag_type == 'dvl'):
            self.upsample_behaviours()
        if (not hasattr(self.upsampled_beh,'beh_simple')) and (self.tag_type == 'dvl'):
            self.test_det_beh_agreement()

        cats = np.unique(self.EthBeh)
        n_cats = len(cats)
        if cols is None:
            # generate distinct colours
            cols = distinctipy.get_colors(n_cats)
        if len(cols) != n_cats:
            raise ValueError(
                f"Number of provided colours ({len(cols)}) does not match the number of detected behaviours ({len(cats)})"
            )
        # generate behaviour-based line collections
        inds, arcs, behavs = self.get_lines_from_string_list(
            string_list=self.EthBeh, line_list=getattr(self.acc, acc_sig)
        )
        # assign relevant colours
        arc_colours = []
        for x in behavs:
            arc_colours.append(list(compress(cols, cats == x))[0])

        fig, ax = plt.subplots(figsize=(6.4, 3.2))
        # set axes limits manually because Collections do not take part in autoscaling
        ax.set_xlim(0, len(self.EthBeh))
        # ax.set_ylim(-6, 6)

        line_collection = LineCollection(arcs, linewidths=1, colors=arc_colours)

        ax.add_collection(line_collection)
        # get order of behaviours and sort for legend
        beh_order = [behavs.index(x) for x in cats]
        # create list of behaviours by order of appearance
        sorted_behaviours = [x for _, x in sorted(zip(beh_order, cats))]
        sorted_cols = [x for _, x in sorted(zip(beh_order, cols))]
        ax.set_ylim(
            np.min(getattr(self.acc, acc_sig)), np.max(getattr(self.acc, acc_sig))
        )
        # generate legend objects
        proxies = [self.make_proxy(x, linewidth=1) for x in sorted_cols]
        # overlay detected foraging if requested
        if plot_vid_forage:
            np.unique(self.upsampled_beh.Behaviour)
            c = np.concatenate(
                (
                    np.where(self.upsampled_beh.Behaviour == "s"),
                    np.where(self.upsampled_beh.Behaviour == "d"),
                ),
                axis=1,
            )
            c.sort(kind="mergesort")
            ax.plot(c.tolist()[0], [1] * c.shape[1], "k*")
            legend_obj = Line2D([], [], color="k", marker="*", linestyle="None")
            proxies.append(legend_obj)
            ax.plot(
                np.where(self.upsampled_beh.beh_simple == "AT")[0],
                [-1] * sum(self.upsampled_beh.beh_simple == "AT"),
                "r*",
            )
            legend_obj = Line2D([], [], color="r", marker="*", linestyle="None")
            proxies.append(legend_obj)
            sorted_behaviours = sorted_behaviours + ["Video foraging", "AT"]
        leg = ax.legend(proxies, sorted_behaviours)
        # increase the width of the legend lines for legelibility
        for line in leg.get_lines():
            line.set_linewidth(4.0)
        ax.set_ylabel(acc_sig)
        ax.title.set_text(f"{self.tagname} estimated behaviours")

        plt.show()

    def test_det_beh_agreement(self, only_forage: bool = False):
        # convert all foraging into 'Forage' and all flight into 'FL'
        self.upsampled_beh["beh_simple"] = self.upsampled_beh.Behaviour
        # simplify foraging
        self.upsampled_beh.loc[
            (self.upsampled_beh.beh_simple == "s")
            | (self.upsampled_beh.beh_simple == "d")
            | (self.upsampled_beh.beh_simple == "CT"),
            "beh_simple",
        ] = "Forage"
        chg_idx, chg_str = self.get_changes_in_string_list(
            self.upsampled_beh.beh_simple
        )
        # define the behaviours found in video
        cats = np.unique(self.upsampled_beh.beh_simple)

        not_AT = np.where(self.upsampled_beh.beh_simple != "AT")[0]

        true_positive_forage_sample_rate = sum(
            (self.EthBeh[not_AT] == "Forage")
            * (self.upsampled_beh.beh_simple[not_AT] == "Forage")
        ) / sum((self.EthBeh[not_AT] == "Forage")).item()

        false_positive_forage_sample_rate = sum(
            (self.EthBeh[not_AT] == "Forage")
            * (self.upsampled_beh.beh_simple[not_AT] != "Forage")
        ) / sum((self.EthBeh[not_AT] == "Forage")).item()

        # how many predicted forages have known foraging in them
        pred_chg_idx, pred_chg_str = self.get_changes_in_string_list(
            self.upsampled_beh.beh_simple
        )

        # which points are consecutive foraging
        pred_forage_changes = np.where(pred_chg_str == "Forage")[0]
        foraging_present = []

        # make copy of ethogram and remove values where bird pecked at tag
        ethocopy = self.EthBeh
        ethocopy[self.upsampled_beh.beh_simple == "AT"] = 0
        for idx in pred_forage_changes:
            if idx == 0:
                ethocopy[0 : (pred_chg_idx[idx] + 1)] == "Forage"
            else:
                foraging_present.append(
                    any(
                        ethocopy[
                            (pred_chg_idx[(idx - 1)] - 1) : (pred_chg_idx[idx] + 1)
                        ]
                        == "Forage"
                    )
                )
                # if any(self.upsampled_beh.beh_simple[(pred_chg_idx[(idx-1)]-1):(pred_chg_idx[idx]+1)] == 'AT'):
                #     foraging_present[-1] = np.nan

        # remove NaNs
        # foraging_present = list(compress(foraging_present,~np.isnan(foraging_present)))
        # non_AT_bouts = len(foraging_present)
        boutTPR = sum(foraging_present) / len(foraging_present)
        boutFPR = sum([not x for x in foraging_present]) / len(foraging_present)

        return (
            true_positive_forage_sample_rate,
            false_positive_forage_sample_rate,
            boutTPR,
            boutFPR,
        )
