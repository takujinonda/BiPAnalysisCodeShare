import pandas as pd
import re


def dtFormat(x):
    """Format incoming datetimes from BiP system.
    Typically datetimes are brought in the format YYYY-mm-dd HH:MM:SS+00:00,
    however, SS can contain decimal values or not. This function returns the
    same datetime in string format with a decimal place added if none is
    originally present and removes '+00:' from the string so that datetimes can
    be converted to datetime format in Pandas.

    Parameters
    ----------
    x
        datetime in string format

    Returns
    -------
    String of x in correct datetime format %Y-%m-%d %H:%M:%S.%f

    """

    # extract all before '+' if present
    if re.findall("^(.*)\+", x):
        out = re.findall("^(.*)\+", x)[0]
    else:
        out = x
    # remove all non-datetime characters
    out = re.sub("[^0-9/ .:-]", "", out)
    if "." not in out:  # add milliseconds if not present (3 d.p.)
        out = out + ".000"
    if out.count(".") > 1:
        out = re.findall("[^.]*.[^.]*", out)[0]  # extract all before second period
    # ensure 3 d.p.
    msLength = len(re.findall("[^.]*.", out)[1])
    if msLength != 3:
        if msLength < 3:
            out = out + ("0" * (3 - msLength))
        else:
            out = out[: -(msLength - 3)]

    return out


def readAxy(
    filename,
    delim="\t",
    cols=None,
    datetimeFormat="%d/%m/%Y %H:%M:%S.%f",
):
    """
    Read in AxyTrek GPS data (txt files) as output by X Manager

    Parameters
    ----------
    filename
        path to AxyTrek txt file
    cols:
        string to denote if acceleration and GPS ('acc') or GPS only ('gps') should be read in
    colnames
        list of string column names to be assigned. Must be of same length as cols. Defaults to ['Date', 'Time', 'lat', 'lon']

    Returns
    -------
    If 'gps' cols passed, single Pandas dataframe of datetime, latitude, and
    longitude.
    If 'acc' cols passed, two Pandas dataframes returned, first an acceleration
    dataframe of formatted datetime (DT), latitude (lat), longitude (lon),
    longitudinal acceleration (X), lateral acceleration (Y), and dorsoventral
    acceleration (Z). Second, a 'gps' formatted dataframe as above.
    If no cols argument given, all data is read in, including pressure
    (pressure), temperature (temp), height about sea level (altitude), and
    ground speed (spd).
    """

    accCols = ["Timestamp", "X", "Y", "Z"]
    gpsCols = ["Timestamp", "location-lat", "location-lon"]

    # read in based on requested columns
    if cols.lower() == "gps":
        try:
            df = pd.read_csv(filename, sep=delim, header=0, usecols=gpsCols)
            df.dropna(inplace=True)
            df.reset_index(inplace=True)
        except ValueError as e:
            if "Usecols do not match columns" in e:
                df = pd.read_csv(filename, sep=",", header=0, usecols=gpsCols)
                df.dropna(inplace=True)
                df.reset_index(inplace=True)
        df.rename(columns={"location-lat": "lat", "location-lon": "lon"}, inplace=True)
    elif cols.lower() == "acc":
        try:
            df = pd.read_csv(
                filename, sep=delim, header=0, usecols=list(set(accCols + gpsCols))
            )
        except ValueError as e:
            if "Usecols do not match columns" in e:
                df = pd.read_csv(
                    filename, sep=",", header=0, usecols=list(set(accCols + gpsCols))
                )
    else:
        df = pd.read_csv(filename, sep=",", header=0)

    df["DT"] = [dtFormat(x) for x in df.Timestamp]  # ensure correct datetime formats
    df["DT"] = pd.to_datetime(df["Timestamp"], format=datetimeFormat)

    # split into acc/gps if acc requested
    if cols == "acc":
        new_acc_cols = [
            "DT",
            "lat",
            "lon",
            "X",
            "Y",
            "Z",
        ]
        new_gps_cols = [
            "DT",
            "lat",
            "lon",
        ]
        accdf = df[new_acc_cols]
        gpsdf = df[new_gps_cols]
        gpsdf.dropna(inplace=True).reset_index(inplace=True)

        return accdf, gpsdf

    else:

        return df


def readBiP(filename, cols=None):
    """
    Read in BiP-formatted data

    Parameters
    ----------
    filename
        string of full path to file.
    cols
        string value indicating if number of columns to read in can be reduced.
        For GPS and acceleration, use 'acc', for GPS only, use 'gps'. Defaults
        to None.

    Returns
    -------
    If 'gps' cols passed, single Pandas dataframe of datetime, latitude, and
    longitude.
    If 'acc' cols passed, two Pandas dataframes returned, first an acceleration
    dataframe of formatted datetime (DT), latitude (lat), longitude (lon),
    longitudinal acceleration (X), lateral acceleration (Y), and dorsoventral
    acceleration (Z). Second, a 'gps' formatted dataframe as above.
    If no cols argument given, all data is read in, including pressure
    (pressure), temperature (temp), height about sea level (altitude), and
    ground speed (spd).
    """

    accCols = [
        "time",
        "acceleration_longitudinal",
        "acceleration_lateral",
        "acceleration_dorso_ventral",
    ]
    gpsCols = [
        "time",
        "latitude",
        "longitude",
    ]

    # read in based on requested columns
    if cols.lower() == "gps":
        df = pd.read_csv(filename, sep=",", header=0, usecols=gpsCols)
        df.dropna(inplace=True)
        df.reset_index(inplace=True)
    elif cols.lower() == "acc":
        df = pd.read_csv(
            filename, sep=",", header=0, usecols=list(set(accCols + gpsCols))
        )
    else:
        df = pd.read_csv(filename, sep=",", header=0)

    # df = pd.read_csv(filename, sep = ",", header = 0, usecols = cols)
    # rename columns for later use
    df.rename(
        columns={
            "time": "DT",
            "latitude": "lat",
            "longitude": "lon",
            "acceleration_longitudinal": "X",
            "acceleration_lateral": "Y",
            "acceleration_dorso_ventral": "Z",
            "pressure": "pressure",
            "temperature": "temp",
            "height_above_mean_sea_level": "altitude",
            "ground_speed": "spd",
        },
        inplace=True,
    )

    df["DT"] = [dtFormat(x) for x in df.DT]  # ensure correct datetime formats
    df["DT"] = pd.to_datetime(df["DT"], format="ISO8601")

    # split into acc/gps if acc requested
    if cols == "acc":
        new_acc_cols = [
            "DT",
            "lat",
            "lon",
            "X",
            "Y",
            "Z",
        ]
        new_gps_cols = [
            "DT",
            "lat",
            "lon",
        ]
        # in case of missing acceleration records, remove observations
        df.dropna(subset=["X","Y","Z"], inplace=True)
        
        accdf = df[new_acc_cols]
        accdf.reset_index(inplace=True)
        gpsdf = df[new_gps_cols]
        gpsdf.dropna(inplace=True)
        gpsdf.reset_index(inplace=True)

        return accdf, gpsdf

    else:

        return df


def readDVL(filename, accStart, fs, vidStart=None, vidOnlyPeriod=False):
    """
    Read Little Leonardo DVL data logger (400M) data. File locations and start-times of acceleration and video recording (dd/mm/yyyy HH:MM:SS) required.
    """
    accSt = pd.to_datetime(accStart)
    vidStart = pd.to_datetime(vidStart)
    # read in acceleration data
    dat = pd.read_table(filename, skiprows=7, sep=",", usecols=[0, 1, 2])

    # remove whitespace from headers
    dat.rename(columns=lambda x: x.strip(), inplace=True)

    # add time series
    dat["DT"] = pd.date_range(accSt, periods=len(dat), freq=f"{int(1/fs*1000)}ms")

    if vidOnlyPeriod:
        # select data within video range
        dat = dat[(dat.DT >= vidStart) & (dat.DT < vidStart + pd.Timedelta(hours=2))]

    return dat.reset_index()
