import platform
import pandas as pd
import glob
import sys
import json
import argparse

from typing import Union
from avimove import birdTag

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

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--dvl-config",
        type=str,
        default="example_DVL_config.json",
        help="File path to config file for DVL files",
    )
    args = parser.parse_args()

    return args

def dvl_from_json(
    filepath: str,
    use_OS: bool = False,
    ):
    """
    Unpack DVL config file and build DVL objects.

    Parameters
    ----------
    filepath: str
        Path to the configuration JSON file.
    use_OS: bool
        If multiple directory paths are given dependent on OS. Keys to path in
        JSON file should conform to platform.system() conventions.
    """
    with open(filepath) as f:
        tag_info = json.load(f)

    if use_OS:
        sys_platform = platform.system()

    # extract the video behaviour classifications
    if use_OS:
        behav_class_file = tag_info["class_folder"][sys_platform]
    else:
        behav_class_file = tag_info["class_folder"]

    if use_OS:
        acc_file = glob.glob(tag_info["parent_folder"][sys_platform] + '*/*acc*.txt')
    else:
        acc_file = glob.glob(tag_info["parent_folder"] + '*/*acc*.txt')

    # extract the tag names, acceleration and video recording start times
    dvls = {}
    for key,value in tag_info["tags"].items():
        # locate the acceleration record
        if use_OS:
            acc_file = glob.glob(tag_info["parent_folder"][sys_platform]+"/"+key+"/*acc*.txt")[0]
        else:
            acc_file = glob.glob(tag_info["parent_folder"]+"/"+key+"/*acc*.txt")[0]
        dvl_config = {
            "tagname": key,
            "accStart": value["acc_start"],
            "vidStart": value["vid_start"],
            "filepath": acc_file,
            "tag_type": "dvl",
            "accfs": 20,
            "long_acc_name": value["acc_names"],
            "gps_fixes_per_minute": None,
        }
        dvls.update(
            {
                key: birdTag(**dvl_config)
            }
        )

    return dvls

dvls=dvl_from_json('avimove/example_DVL_config.json',use_OS=True)

for x in dvls.values():
    x.readin(vidOnlyPeriod=True)
    x.readBeh(behav_class_loc)
    x.upsample_behaviours()
    x.flight_est(numPoints=10,removeErr=True)
    x.flapping()
    x.calculate_thresholds()

if all(b == 'dvl' for b in [x.tag_type for x in dvls.values()]):
    # calculate pitch differences for all tags
    medianPitchChanges = []
    for x in dvls.values():
        medianPitchChanges.append(x.calculate_thresholds(threshold_scale = 1.5))
    from statistics import median
    toEx = median(medianPitchChanges)


dvls['17008'].plot_acc_behaviours('pitch')

import numpy as np
# define an empty ethogram
dvls['17008'].EthBeh = np.array(["Unknown" for _ in range(len(dvls['17008'].acc))])
dvls['17008'].EthBeh[np.where(dvls['17008'].acc.ODmn < 0.2)] = "Rest"
DiveUp = median([np.mean(dvls['17008'].acc.pitch[x]) for x in dvls['17008'].flInds]) + 2 * median(
    [np.var(dvls['17008'].acc.pitch[x]) for x in dvls['17008'].flInds]
)
DiveDown = median([np.min(dvls['17008'].acc.pitch[x]) for x in dvls['17008'].flInds]) - 30
dvls['17008'].EthBeh[np.where(dvls['17008'].flap_bouts == 1)] = "FL"
# find pitch changes
pitpeaks, pittroughs, _ = accFn.peak_trough(dvls['17008'].acc.pitch)
PitUp = np.array(dvls['17008'].acc.pitch[pitpeaks]) - np.array(
    dvls['17008'].acc.pitch[pittroughs]
)
PitDown = np.array(dvls['17008'].acc.pitch[pittroughs]) - np.array(
    dvls['17008'].acc.pitch[pitpeaks]
)
PitUpL = list(compress(pitpeaks, abs(PitDown) > toEx))
PitDownL = list(compress(pittroughs, PitUp > toEx))
PitLarge = sorted(PitUpL + PitUpL)


sum((test.EthBeh[not_AT] == "Forage") * (test.upsampled_beh.beh_simple[not_AT] == "Forage")) / sum((test.EthBeh[not_AT] == "Forage")).item()