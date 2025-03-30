# Script to assess accuracy of behaviour detection method on DVL tag data
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

def run_behaviour_detection(tag_data: birdTag):
    """
    Run through the behaviour detection system of a birdTag object.

    Returns the true and false positive rates on sample and bout levels (the
    number of samples/bouts correctly/falsely predicted, respectively).

    Parameters
    ----------
    tag_data
        birdTag object of tags to be analysed.

    Returns
    -------
    accuracies
        Dictionary of TP/FP rates for each tag with tagname as key.
    """
    for x in tag_data:
        print(f'Reading in {x.tagname} acc data...')
        x.readin(vidOnlyPeriod=True)
        print('Finished reading in acc.')
        print(f'Reading in {x.tagname} behaviours...')
        x.readBeh(behav_class_loc)
        print('Finished reading in behaviours.')
        print(f'Upsampling behaviour data...')
        x.upsample_behaviours()
        print('Finished regular behaviours.')
        print('Calculating flight estimated periods...')
        x.flight_est(numPoints=10,removeErr=True)
        print('Finished flight estimates.')
        print('Calculating flapping periods')
        x.flapping()
        print('calculate flight pitch changes')
        x.calculate_thresholds()
        print('Thresholds calculated')

    if all(b == 'dvl' for b in [x.tag_type for x in tag_data]):
        # calculate pitch differences for all tags
        medianPitchChanges = []
        for x in tag_data:
            medianPitchChanges.append(x.calculate_thresholds(threshold_scale = 1.5))
        from statistics import median
        toEx = median(medianPitchChanges)

    accuracies = {}
    for tag in tag_data:
        tag.beh_detect(toEx=toEx)
        sTPR,sFPR,bTPR,bFPR = tag.test_det_beh_agreement()
        accuracies.update(
            {
                tag.tagname:[sTPR,sFPR,bTPR,bFPR]
            }
        )

    # convert to pd DataFrame
    accuracies = pd.DataFrame(data=accuracies,index=['sTRP','sFPR','bTPR','bFPR']).transpose()

    return accuracies

def main():
    args = parse_args()
    dvls=dvl_from_json(args.dvl_config,use_OS=True)
    accuracies = run_behaviour_detection(list(dvls.values()))
    print(accuracies)

if __name__=="__main__":
    main()