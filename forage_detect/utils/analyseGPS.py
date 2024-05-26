def gps_distanceSingle(longitudes, latitudes, latVal, lonVal):                                                                       
    # taken from https://www.tjansson.dk/2021/03/vectorized-gps-distance-speed-calculation-for-pandas/, thanks to Thomas Jansson for use of this function
    """                                                                                                                 
    Calculates distances between an array of GPS points and a single position via a vectorized haversine calculation the great circle distance between two arrays of points on the earth (specified in decimal degrees). All args must be of equal length.                                         
 
    Args:                                                                                                               
        longitudes: pandas series of longitudes                                                                         
        latitudes:  pandas series of latitudes                               
        latVal:     latitude value of single position                                           
        lonVal:     longitude value of single position
 
    Returns:                                                                                                            
        Distance is returned an array in meters.                                                                             
 
    Example:                                                                                                            
        >>> df['gpsDist'] = gps_distance(df.longitude, df.latitude, 39.400,141.998)
    """                                                                  

    return [distance.distance((x,y),(latVal,lonVal)).km for x,y in zip(latitudes,longitudes)]

def gps_speed(longitudes, latitudes, timestamps):                                                                       
    # taken from https://www.tjansson.dk/2021/03/vectorized-gps-distance-speed-calculation-for-pandas/, thanks to Thomas Jansson for use of this function
    """                                                                                                                 
    Calculates the instantaneous speed from the GPS positions and timestamps. The distances between the points          
    are calculated using a vectorized haversine calculation the great circle distance between two arrays of points on   
    the earth (specified in decimal degrees). All args must be of equal length.                                         
 
    Args:                                                                                                               
        longitudes: pandas series of longitudes.                                                                         
        latitudes:  pandas series of latitudes.                                                                          
        timestamps: pandas series of timestamps.                                                                         
 
    Returns:                                                                                                            
        Speed is returned an array in m/s.                                                                             
 
    Example:                                                                                                            
        >>> df['gpsSpeed'] = gps_speed(df.longitude, df.latitude, df.recordedAt).
    """
 
    lon1 = longitudes[:-1].reset_index(drop = True)                                                                                       
    lat1 = latitudes[:-1].reset_index(drop = True)                                                                                 
    lon2 = longitudes[1:].reset_index(drop = True)                                                                                       
    lat2 = latitudes[1:].reset_index(drop = True)                                                                                 
 
    # Vectorized haversine calculation                                                                                  
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])                                                  
    a = np.sin((lat2 - lat1) / 2.0)**2 + (np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2.0)**2)                 
    dist = np.array(6371 * 2 * np.arcsin(np.sqrt(a)) * 1000)
    time_array = (timestamps.diff().dt.seconds ).values[1:]                                                       
 
    # Calculate the speed                                                                                               
    time_array[time_array == 0] = np.nan  # To avoid division by zero                                                   
    speed = np.array(dist / time_array)

    # Make the arrays as long as the input arrays                                                                        
    speed = np.insert(speed, 0, np.nan, axis=0)
    dist = np.insert(dist, 0, np.nan, axis=0)                                                               
    return dist,speed

def distSpeed(lat,lon,DT,threshold=None):
    """
    Calculates distances (metres) and speed (m/s) of GPS positions and datetime information. Uses a speed threshold (default None) to define erroneous GPS positions. GPS and datetime lengths must agree.

    Args:

        lat:    array of latitude points in decimal degrees.
        lon:    array of longitude points in decimal degrees.
        DT:     datetime array in Pandas datetime format.

    Returns:
        Float arrays of distance and speed values.
    """

    dist,speed = gps_speed(lat,lon,DT)
    if threshold is not None:

        while np.nanmax(speed) > threshold:

            # remove erroneous GPS values
            lon[speed > threshold] = np.nan
            lat[speed > threshold] = np.nan

            # recalculate speed
            dist,speed = gps_speed(lat,lon,DT)
    
    return dist,speed

def removeNear(data, captureSite, distThreshold = 1.5):
    """
    Remove acceleration and GPS data captured within a set distance of the capture site.

    Args:
        data:           pandas dataframe of tag data as formatted by readBIP. Expected columns are DT (datetime, formatted), lat (latitude in decimal degrees), and lon (longitude in decimal degrees).
        captureSite:    capture site latitude and longitude in decimal degrees.
        distThreshold:  desired distance threshold from capture site in km. Defaults to 1.5.

    Returns:
        Pandas dataframe of similar format as data with all points located within distThreshold of captureSite removed.
    """
    # highest data sample rate (Hz)
    fs = int(np.timedelta64(1,'s') / stats.mode(np.diff(data.DT),keepdims=False)[0])

    # GPS sample rate (fixes per minute)
    gpsFS = 60/int(np.timedelta64(1,'m') / stats.mode(np.diff(data.dropna()['DT']),keepdims=False)[0])

    # find distances from capture site
    capSiteDist = np.array(gps_distanceSingle(data.lon.dropna(),data.lat.dropna(),captureSite[0],captureSite[1]))

    inds = data.lat.dropna().index[np.where(capSiteDist < distThreshold)[0]] # find indeces to be removed from original dataset

    # if any difference in GPS and acc sampling rates means we must extend GPS indeces
    samplerDiff = int(np.floor(fs * gpsFS / 2)) # number of fs samples either side of GPS position

    # create ranges to be removed
    rem = []
    for x in inds:
        rem.extend(np.arange(x - samplerDiff, x + samplerDiff))
    rem = np.array(rem)

    return data.drop(rem)