"""Wrapper script for birdTag functionality.
Use arguments to define which functions are desired, the format of output, and
the data to analyse.
"""

import argparse
from BirdTag import BirdTag
from numpy import seterr


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
        "-tname",
        "--tag-name",
        type=str,
        default=None,
        help="Name of tag (used in generating output file names). \
            Defaults to the base file name.",
    )
    # parser.add_argument(
    #     "-f",
    #     "--forage",
    #     choices=[True,False],
    #     default=False,
    #     help="Estimate foraging behaviours from acceleration records.",
    # )
    parser.add_argument(
        "-fg",
        "--flap-glide",
        action="store_true",
        help="Derive flapping/gliding from acceleration records and save to file?",
    )
    parser.add_argument(
        "-w",
        "--wind",
        action="store_true",
        help="Estimate wind speeds from GPS records?",
    )
    parser.add_argument(
        "-rem-loc",
        "--remove-location",
        nargs=2,
        default=None,
        help="Lat lon of specified location to remove nearby values \
            from (decimal degrees).",
    )
    parser.add_argument(
        "-rem-dist",
        "--remove-distance",
        type=float,
        default=1.5,
        help="Buffer from specified location to remove data (km). Default 1.5.",
    )
    # parser.add_argument(
    #     "-n-fl",
    #     "--n-flight-points",
    #     type=int,
    #     default=24,
    #     help="Number of minutes per day to estimate flight from \
    #         frequency analysis.",
    # )
    parser.add_argument(
        "-spthresh",
        "--speed-threshold",
        type=float,
        default=24,
        help="Speed threshold for GPS-based speed calculations (m/s).\
            Default 24.",
    )
    parser.add_argument(
        "-sf",
        "--saveloc",
        type=str,
        default=None,
        help="Path for output(s). Meaningful file names are automatically generated.\
            If no path provided, outputs saved in current working directory.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Verbosity.",
    )
    args = parser.parse_args()
    return args


def main():
    """Wrapper for full BirdTag functionality.
    Based on CLI arguments, various analyses will be performed with desired
    outputs written to file.
    """
    args = parse_arguments()

    # init birdTag object
    tag = BirdTag(
        filepath=args.data_url,
        tag_type="bip",
        tagname=args.tag_name,
        analysis_outloc=args.saveloc,
        verbose=args.verbose,
    )

    # read in the relevant data
    if (args.wind) * (not args.flap_glide):
        tag.readin(wind=True)
    else:
        tag.readin(wind=False)

    # perform the wind analysis if requested
    if args.wind:
        if args.verbose:
            print("Running wind analysis...")
        seterr(invalid="ignore", over="ignore")  # temporary error suppression
        tag.estimate_wind()

    # calculate flap/glide if requested
    if args.flap_glide:
        # remove data near position if requested
        if args.remove_location:
            if args.verbose:
                print(f"Removing points {args.rem_dist} km from {args.remove_location}...")
            tag.remove_near(args.remove_location, args.remove_distance)
        if args.verbose:
            print("Running flap/glide analysis...")
        tag.save_flap_glide()


if __name__ == "__main__":
    main()
