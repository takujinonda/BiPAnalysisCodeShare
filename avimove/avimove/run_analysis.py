"""Wrapper script for birdTag functionality.
Use arguments to define which functions are desired, the format of output, and
the data to analyse.
"""
import argparse
import birdtag
import json

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-url",
        "--data-url",
        type=str,
        required=True,
        help="URL of the data to be analysed.",
    )
    parser.add_argument(
        "-conf",
        "--config",
        type=str,
        required=True,
        help="Path to tag config file.",
    )
    parser.add_argument(
        "-f",
        "--forage",
        choices=[True,False],
        default=False,
        help="Estimate foraging behaviours from acceleration records.",
    )
    parser.add_argument(
        "-rem-near",
        "--remove-near",
        choices=[True,False],
        default=False,
        help="Remove values near specified location.",
    )
    parser.add_argument(
        "-rem-loc",
        "--remove-location",
        nargs=2,
        default=None,
        help="Lat lon of specified location to remove nearby values from (decimal degrees).",
    )
    parser.add_argument(
        "-n-fl",
        "--n-flight-points",
        type=int,
        default=24,
        help="Number of minute periods to estimate flight from frequency analysis.",
    )
    parser.add_argument(
        "-sf",
        "--savefile",
        type=str,
        default="",
        help="Path and filename pattern for output(s). Meaningful file names \
            are automatically suffixed to pattern. If no path provided, \
            outputs saved in local directory.",
    )
    args = parser.parse_args()
    return args

def main():
    """Wrapper for full birdTag functionality.
    Based on CLI arguments, various analyses will be performed with desired
    outputs written to file.
    """
    args = parse_arguments()
    
    with open(args.config,'r') as f:
        config = json.load(f)

    # init birdTag object
    tag = birdTag(
        filepath = args.data_url,
        **config
    )
    
    # read in the relevant data
    tag.readin()
    # remove data near position if requested
    if args.rem_near:
        tag.remove_near(args.remove_location)
    # perform the relevant analysis as requested
    if args.forage:
        tag.estimate_flight(
            num_points = args.n_flight_points
            )
        tag.flappign
        toEx = tag.calculate_thresholds()
        tag.beh_detect(toEx)
    
    if 