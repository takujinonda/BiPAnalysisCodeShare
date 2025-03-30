from wind import derive_wind as windFn
from read import read_in as read_in
from foraging import acceleration_analysis as accFn

import matplotlib.pyplot as plt
import distinctipy
import numpy as np
import pandas as pd

# suppress warnings on chained assignment (not relevant here)
pd.options.mode.chained_assignment = None
import datetime as dt
from typing import Union

from mpl_toolkits.basemap import Basemap
from statistics import median
from typing import Union
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
from itertools import compress
from pathlib import Path
import os


class BirdTag:
    """Class for avian biologging tags.
    Current support includes AxyTrek with a look toward incorporating NinjaScan
    data. Requires file path(s) and string descriptor of tag ('Axy').
    Methods include reading of data from raw CSV/txt or BiP database. Data loads
    can include GPS or not (where present).
    Vectorised haversine distance and speed calculation from GPS. 
    Dynamic/static acceleration calculation via equiripple lowpass filter. 
    Pitch calculation.
    Data removal from proximity to location.
    """

    def __init__(
        self,
        filepath: str,
        tag_type: str,
        tagname: str = None,
        accfs: Union[int, None] = None,
        acc_name_format: Union[str, None] = None,
        analysis_outloc: Path = None,
        accStart: Union[np.datetime64, None] = None,
        vidStart: Union[np.datetime64, None] = None,
        verbose: bool = False,
        *args,
        **kwargs,
    ):
        """
        Parameters
        ----------
        filepath
            Path to signal data
        tag_type
            Type of tag. Currently choice of 'BiP' and 'AxyTrek' (case
            insensitive)
        tagname
            Tag identifier. If none provided, the file base name is used.
        accfs
            Acceleration sampling frequency (Hz). If not provided, will be
            derived from acceleration data if read in.
        acc_name_format
            list of names for acceleration signals. Must conform to the
            following order - longitudinal (along body axis), dorsoventral
            (vertical axis), lateral
        analysis_outloc
            Path to desired save location of analysis results.
        accStart
            DateTime start of acceleration signal recording
        vidStart
            DateTime start of video recording
        """
        self.filepath = filepath
        self.tag_type = tag_type.casefold()
        if tagname is None:
            tagname = os.path.basename(self.filepath)
        self.tagname = tagname
        self.accfs = accfs
        self.acc_name_format = acc_name_format
        self.accStart = pd.to_datetime(accStart, format="%d/%m/%Y %H:%M:%S")
        if vidStart is not None:
            vidStart = pd.to_datetime(vidStart, format="%d/%m/%Y %H:%M:%S")
        self.vidStart = vidStart
        if analysis_outloc is None:
            analysis_outloc = os.getcwd()
        self.analysis_outloc = analysis_outloc
        self.verbose = verbose

        # define some filenames for analysis output
        self.wind_out_filename = self.tagname + "_estimated_wind.csv"

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

    def readin(self, wind: bool, vidOnlyPeriod: bool = True) -> None:
        if self.verbose:
            print("Reading file...")
        if (self.tag_type == "axy") * wind:
            self.gps = read_in.readAxy(self.filepath, cols="gps")
        elif self.tag_type == "axy":
            self.acc, self.gps = read_in.readAxy(self.filepath, cols="acc")
        if (self.tag_type == "bip") * wind:
            self.gps = read_in.readBiP(self.filepath, cols="gps")
        elif self.tag_type == "bip":
            self.acc, self.gps = read_in.readBiP(self.filepath, cols="acc")
        # derive acceleration sampling frequency if not provided
        if (self.accfs is None) & (hasattr(self, "acc")):
            self.accfs = np.timedelta64(1, "s") / np.mean(np.diff(self.acc.DT))

    def remove_near(self, home_site, dist_threshold=1.5) -> None:
        """Identify and remove data near a specified location.
        Both acceleration and GPS data will be removed for
        positions/accelerometer data found within threshold distance of
        home_site.

        Parameters
        ----------
        home_site
            Location to find nearby data (lat, lon - decimal degrees).
        dist_threshold
            Distance within which to identify data points.
        """
        acc_inds = windFn.removeNear(
            acc_data=self.acc,
            gps_data=self.gps,
            captureSite=home_site,
            distThreshold=dist_threshold,
        )

        # remove acceleration data
        self.acc.drop(acc_inds, inplace=True)
        self.acc.reset_index(inplace=True)
        # repeat for GPS data
        self.gps.drop(set(acc_inds).intersection(set(self.gps.index)), inplace=True)
        self.gps.reset_index(
            inplace=True, drop=True
        )  # gps data already has an index reference

    def dist_speed(self, threshold=None) -> None:
        """Calculate distance and speed with optional maximum speed threshold.
        Speed and distance added as attributes.

        Parameters
        ----------
        threshold
            Optional maximum speed threshold (m/s). If exceeded, suspect GPS
            positions are removed and speed recalculated until no more threshold
            breaches remain.
        """
        self.dist, self.speed = windFn.distSpeed(
            lat=self.gps.lat, lon=self.gps.lon, DT=self.gps.DT, threshold=threshold
        )

    def estimate_wind(self) -> None:
        """Wrapper function for wind estimation method."""
        wind_df = windFn.windEstimation2(self.gps)
        wind_out_filename = (
            self.analysis_outloc + "/" + self.tagname + "_wind_estimation.csv"
        )
        if self.verbose:
            print(f"Saving wind results to {wind_out_filename}")

        wind_df.to_csv(Path(wind_out_filename), date_format="%Y-%m-%d %H:%M:%S.%f")

    @staticmethod
    def round_seconds(obj: dt.datetime) -> dt.datetime:
        """
        Round a datetime to the nearest second.
        """
        if obj.microsecond >= 500_000:
            obj += dt.timedelta(seconds=1)
        return obj.replace(microsecond=0)

    def accFeatures(self, passband: float = 1.5, stopband: float = 3.0) -> None:
        """Derive acceleration features (dynamic/static).
        Generates dynamic and static acceleration values and outputs dynamic
        surge and dorsoventral accelerations, ODBA, pitch angles, and their
        moving means.

        Parameters
        ----------
        passband
            Passband for low pass filtering for static acceleration (Hz). Default 1.5.
        stopband
            Stopband for low pass filtering for static acceleration (Hz). Default 3.
        """
        # check if data read in
        if not hasattr(self, "acc"):
            print("No acceleration data present")
        else:
            out = accFn.accFeatures(
                self.acc, self.acc_name_format, passband, stopband, self.accfs
            )
            out.index = self.acc.index
            self.acc = pd.concat([self.acc, out], axis=1)

    def roll_sum(self, min_freq: float = 3.0, max_freq: float = 5.0) -> None:
        """Calculate rolling sum of dorsoventral acceleration frequency content.
        Frequency spectral intensity summed between min and max frequencies and
        summed over minute values. Used in generating flight estimates, so
        minimum and maximum frequencies should reflect expected range of
        species' flapping.

        Parameters
        ----------
        min_freq
            Minimum frequency of summation (Hz). Default 3.0.
        max_freq
            Maximum frequency of summation (Hz). Default 5.0.
        """
        self.rolling_freq_sum = accFn.rollingSpecSum(
            sig=self.acc.Z,
            min_freq=min_freq,
            max_freq=max_freq,
            fs=self.accfs,
            dur=60,
            inclusive=False,
        )

    def flight_est(self, num_points: int) -> None:
        """Estimate flight from frequency content of dorsoventral acceleration.
        Using summed spectral intensity of expected flapping frequencies, derive
        estimates of flight periods separated by a minimum time gap.

        Parameters
        ----------
        num_points
            Number of flapping periods to estimate per full day of recording.
        """
        # find range of full days in dataset
        days_of_accel_data = round(
            (self.acc.DT.iloc[0] - self.acc.DT.iloc[0]) / np.timedelta64(1, "D")
        )
        num_points = num_points * days_of_accel_data

        # check if roll_sum calculation has been performed
        if not hasattr(self, "rolling_freq_sum"):
            if self.verbose:
                print("Performing spectral rolling sum")
            self.roll_sum()
            
        self.flInds, self.flight = accFn.flightestimate(
            signal=self.acc.Z,
            roll_sum=self.rolling_freq_sum,
            fs=self.accfs,
            dt=self.acc.DT,
            num_points=num_points,
        )

    def flight_thresholds(self) -> None:
        if (not hasattr(self, "flInds")) | (not hasattr(self, "acc")):
            raise AttributeError("Flight estimate or acceleration data required")
        else:
            self.acc = pd.concat(
                [self.acc, accFn.flight_est_thresholds(self.acc, self.flInds)], axis=1
            )

    @staticmethod
    def flatten(xss):
        return [x for xs in xss for x in xs]

    def flapping(self) -> None:
        """Estimate flapping and gliding behaviour.
        Calculate mask of flapping behaviour and 'bouts', grouped by `flap_freq`
        and `bout_gap` seconds, respectively. This function requires a flight
        mask to be present.
        Glides are assigned during flap bouts but when flapping not assigned.
        """

        (
            self.flap,
            self.flap_bouts,
            self.flap_start,
            self.flap_end,
            self.flap_peaks,
            self.flap_troughs,
        ) = accFn.flap(sig=self.acc.Z, fs=self.accfs, bout_gap=10, flap_freq=4)
        # define a glide mask too
        self.glide = ((self.flap == 0) & (self.flap_bouts == 1)).astype(int)

    def pitchPT(self) -> None:
        """
        Calculate pitch thresholds from estimated flight periods
        """
        if not hasattr(self, "flInds"):
            raise AttributeError("Add estimated flight periods")

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
                self.roll_sum()
            if not hasattr(self, "flight"):
                print("Flight not calculated\nRunning now...")
                self.flight_est()

        if (toEx is None) & (self.tag_type == "axy"):
            raise ValueError("Pitch threshold not given")

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

    def save_flap_glide(
        self,
        sp_threshold: float = 80,
    ):
        """Save flapping/gliding and speed information to file.
        Wrapper function to put together analysis of dorsoventral acceleration
        to derive flapping and glide data with speeds.

        Parameters
        ----------
        sp_threshold
            Speed threshold (m/s) to remove erroneous values. Default 80.
        """
        self.flapping()
        self.dist_speed(sp_threshold)

        # create GPS series of same length as acceleration
        if "index" in self.acc.columns:
            acc_idx = self.acc["index"]
            gps_idx = self.gps["index"]
        else:
            acc_idx = self.acc.index
            gps_idx = self.gps.index
        dummy_series = pd.Series(np.nan, acc_idx)
        gps_data = {
            "lat": dummy_series.copy(),
            "lon": dummy_series.copy(),
            "speed": dummy_series.copy(),
        }
        for key in gps_data.keys():
            if key == "speed":
                gps_data[key][gps_idx] = self.speed
            else:
                gps_data[key][gps_idx] = self.gps[key]

        # combine information into single database
        flap_glide_out = pd.concat(
            [
                self.acc.DT,
                pd.Series(self.flap),
                pd.Series(self.flap_bouts),
                pd.Series(self.glide),
            ],
            axis=1,
            ignore_index=True,
        )
        flap_glide_out = pd.concat(
            [
                flap_glide_out,
                pd.DataFrame(data=gps_data).set_index(flap_glide_out.index),
            ],
            axis=1,
        )

        # create filename and save
        flap_glide_filename = (
            self.analysis_outloc + "/" + self.tagname + "_flap_glide_speed.csv"
        )
        if self.verbose:
            print(f"Saving flap/glide results to {flap_glide_filename}")

        flap_glide_out.to_csv(
            Path(flap_glide_filename), date_format="%Y-%m-%d %H:%M:%S.%f"
        )

    def calculate_thresholds(
        self,
        threshold_scale: float = 1.5,
        ) -> None:
        """
        Calculate all the required threshold for behaviour estimation.
        """
        # check if acceleration characteristics are present
        if not "pitch" in self.acc:
            self.accFeatures()
        # check if flight estimate periods are present
        if not hasattr(self, "flInds"):
            print("Estimating flight periods: number of estimated flight minutes requested \n")
            num_points = input("Minutes of flight to estimate per day: ")
            self.flight_est(num_points=num_points)
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

    def plot_acc_behaviours(self, acc_sig, cols=None, plot_vid_forage: bool = False):

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
            string_list=self.EthBeh, line_list=self.acc[acc_sig]
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
        ax.set_ylim(np.min(self.acc[acc_sig]), np.max(self.acc[acc_sig]))
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

    def plot_forage_locations(
        self,
        projection: str = "gall",
        buffer_lat: int = 5,
        buffer_lon: int = 5,
    ):
        """Plot predicted foraging locations.
        Ethogram is required to run this. A global map coastline will be plotted
        and the plot will show the max/min lat lon values + some optional buffer
        (default 5).

        Parameters
        ----------
        projection
            Basemap projection string.
        buffer_lat
            Buffer for latitude plotting in decimal degrees. Default 5.
        buffer_lon
            Buffer for longitude plitting in decimal degrees. Default 5.
        """
        m = Basemap(
            projection=projection,
        )

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

        true_positive_forage_sample_rate = (
            sum(
                (self.EthBeh[not_AT] == "Forage")
                * (self.upsampled_beh.beh_simple[not_AT] == "Forage")
            )
            / sum((self.EthBeh[not_AT] == "Forage")).item()
        )

        false_positive_forage_sample_rate = (
            sum(
                (self.EthBeh[not_AT] == "Forage")
                * (self.upsampled_beh.beh_simple[not_AT] != "Forage")
            )
            / sum((self.EthBeh[not_AT] == "Forage")).item()
        )

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
