from avimove import birdTag
from avimove.forage_detect.src.utils.loadIn import readAxy, dtFormat
import json


# json_file = "avimove/axytrek_config.json"

# with open(json_file,'r') as f:
#     axy_details = json.load(f)


# file_loc = '/Users/aran.garrod/Library/CloudStorage/GoogleDrive-aran.garrod@googlemail.com/My Drive/Data/AxyTrek/2017-09/2017-9_S1.csv'
file_loc = "/Users/aran.garrod/Downloads/ed9ff124-588d-4053-80a6-a326500ce49b_std.csv"
axy_file_config = {
    "filepath": file_loc,
    "tag_type": "bip",
    "tagname": "2017-9",
    "accfs": 25,
    "acc_name_format":["X","Z","Y"],
    "gps_fixes_per_minute": 1,
}

test=birdTag(**axy_file_config)
test.readin()
test.remove_near([39.400,141.998])
test.flight_est(num_points=24*4)
toEx = test.calculate_thresholds()
test.flapping()
test.dist_speed(80)
test.beh_detect(toEx)

a = ['DT',
'lat',
'lon',
'X',
'Y',
'Z',]

# create GPS series of same length as acceleration
if 'index' in test.acc.columns:
    acc_idx = test.acc['index']
    gps_idx = test.gps['index']
else:
    acc_idx = test.acc.index
    gps_idx = test.gps.index
dummy_series = pd.Series(np.nan, acc_idx)
gps_data = {
    'lat': dummy_series.copy(),
    'lon': dummy_series.copy(),
    'speed': dummy_series.copy()}
for key in gps_data.keys():
    if key == 'speed':
        gps_data[key][gps_idx] = test.speed
    else:
        gps_data[key][gps_idx] = test.gps[key]

test_out = pd.concat(
    [
        test.acc.DT,
        pd.Series((test.EthBeh == "Rest") & (test.flight == 0)),
        pd.Series(test.flap),
        pd.Series(test.flap_bouts),
        pd.Series(test.glide),
        pd.Series(test.gps.lat),
        pd.Series(test.gps.lon),
        pd.Series(test.EthBeh),
        ],
    axis=1
    )
test_out.iloc[:,'Speed'] = []
test_out.rename(columns={
    0:'on_water',
    1:'flapping',
    2:'flapping_bouts',
    3:'gliding',
    4:'lat',
    5:'lon',
    6:'Ethogram'},
                inplace=True
                )

test_out.to_csv('example_output.csv',date_format="%Y-%m-%d %H/%M/%S.%f")

import matplotlib.pyplot as plt
plt.plot(test.acc.DT,test.acc.Z,)
for segment in test.flInds:
    plt.plot(
        test.acc.DT[segment],
        test.acc.Z[segment],
        'r',
        )
plt.show()

test.gps



acc_data = test.acc
gps_data = test.gps
captureSite = [39.400,141.998]
from scipy import stats
fs = int(np.timedelta64(1,'s') / stats.mode(np.diff(acc_data.DT),keepdims=False)[0])

# GPS sample rate (fixes per minute)
gpsFS = 60/int(np.timedelta64(1,'m') / stats.mode(np.diff(gps_data['DT']),keepdims=False)[0])

# find distances from capture site
capSiteDist = np.array(gpsFn.gps_distanceSingle(gps_data.lon,gps_data.lat,captureSite[0],captureSite[1]))

inds = gps_data.lat.index[np.where(capSiteDist < distThreshold)[0]] # find indices to be removed from original dataset

# if any difference in GPS and acc sampling rates means we must extend GPS indices
samplerDiff = int(np.floor(fs * gpsFS / 2)) # number of fs samples either side of GPS position

# create ranges to be removed
rem = []
for x in inds:
    rem.extend(np.arange(x - samplerDiff, x + samplerDiff))
rem = np.array(rem)

rem = test.id_near([39.400,141.998])
test.acc.drop(rem,inplace=True).reset_index(inplace=True)
test.gps = test.gps.drop(rem).reset_index()


import pandas as pd
gpsCols = ['Timestamp','location-lat','location-lon']
filename=test.filepath
delim='\t'
try:
    df = pd.read_csv(filename, sep = delim, header = 0, usecols = gpsCols)
    df.dropna(inplace=True)
except ValueError as e:
    if "Usecols do not match columns" in e:
        df = pd.read_csv(filename, sep = ",", header = 0, usecols = gpsCols).dropna(inplace=True)
df.rename(columns={'location-lat':'lat','location-lon':'lon'},inplace=True)
df = pd.read_csv(test.filepath,sep='\t',header=0,usecols=['Timestamp','location-lat','location-lon'])
df.dropna(inplace=True)


blah = test.acc[['DT','Z']]
blah['Speed'] = pd.Series(dtype='float')
blah.loc[test.gps['index'],'Speed'] = test.speed.tolist()


blah.loc[:4,'Speed'] = [1,2,3,4,5]




plt.plot(
    test.acc.DT,
    test.acc.DZ
)
plt.plot(
    test.acc.DT[(test.flap == 1) & (test.EthBeh == "FL")],
    test.acc.DZ[(test.flap == 1) & (test.EthBeh == "FL")],
    'r+',
)
plt.plot(
    test.acc.DT[(test.flap_bouts == 1) & (test.EthBeh == "FL")],
    [0] * sum((test.flap_bouts == 1) & (test.EthBeh == "FL")),
    'g+'
)
plt.show()


for idx in test.flInds:
    plt.scatter(test.acc.DT[test.acc.index[idx]],[1]*len(idx),s=5)
plt.plot(test.gps.DT,test.speed)
plt.show()

fs=test.accfs
sig=test.rolling_freq_sum
window=60
min_gap=5
num_points=20


b = np.zeros(shape=(num_points,window*fs)).astype(int)

out = np.zeros(shape=(num_points,window*fs)).astype(int)
sigad = sig.copy()
for idx in range(out.shape[0]):
    out[idx] = np.arange(np.argmax(sigad)*2 - (round(fs*window/2)),
            (np.argmax(sigad)*2 + round(fs*window/2)))
    
    b[idx] = np.arange(np.argmax(sigad) - (round(fs*window/2)),
                  (np.argmax(sigad) + round(fs*window/2)))
    
    # reduce magnitude of this period and 5 surrounding minutes
    sigad[
        np.arange(np.argmax(sigad) - (round(fs*window/2)),
                  (np.argmax(sigad) + round(fs*window/2)))
    ] = np.min(sigad)


b = np.arange(np.argmax(sigad) - (round(fs*window/2)),
          (np.argmax(sigad) + round(fs*window/2)))

fig,ax = plt.subplots()
ax.plot(test.rolling_freq_sum.index,test.rolling_freq_sum)
for segment in b:
    ax.plot(segment,test.rolling_freq_sum[segment],c='r')
plt.show()

plt.plot(test.acc.DT,test.acc.Z)
plt.plot(test.acc.DT[out[0]],test.acc.Z[out[0]],'r')
plt.show()

plt.plot(test.rolling_freq_sum)
plt.show()

glides = (test.EthBeh == "FL") & (test.flap == 0) & (test.flap_bouts == 1)

import numpy as np
import pandas as pd
test.gps.lon != np.nan

import avimove.forage_detect.src.utils.analyseAcc as accFn

accFn.flight_pitch_changes(test.acc.pitch,test.flInds,'max')

for flidx in test.flInds:
    pk,tr,_=accFn.peak_trough(test.acc.pitch[flidx])
    pk=flidx[pk]
    tr=flidx[tr]

    np.min(np.array(test.acc.pitch[pk]) - np.array(test.acc.pitch[tr]))


    np.max(np.array(test.acc.pitch[pk]) - np.array(test.acc.pitch[tr]))

plt.plot(test.acc.pitch[test.flInds[0]])
plt.plot(
    pk,
    test.acc.pitch[test.flInds[0]][pk],
    'r+',
    )
plt.plot(
    tr,
    test.acc.pitch[test.flInds[0]][tr],
    'g+',
)
plt.show()

b = find_peaks(test.acc.pitch[test.flInds[0]])[0]

from statistics import median
median(np.array(test.acc.pitch[tr]))

from scipy.signal import find_peaks
import numpy as np
import matplotlib.pyplot as plt

plt.plot(
    test.acc.DT,
    test.acc.Z
)
plt.plot(
    test.acc.DT[test.flap_peaks],
    test.acc.Z[test.flap_peaks],
    'r+',
)
plt.plot(
    test.acc.DT[test.flap_troughs],
    test.acc.Z[test.flap_troughs],
    'g+',
)
plt.plot(
    test.acc.DT[test.glide==1],
    test.acc.Z[test.glide == 1],
    'k+',
)
plt.show()

sig = test.acc.Z
peaks,_ = find_peaks(sig)
troughs,_ = find_peaks(-sig)

while len(peaks) != len(troughs):
    if peaks[0] > troughs[0]:
        troughs = np.delete(troughs,0)
    if peaks[-1] > troughs[-1]:
        peaks = np.delete(peaks,-1)

# create alternative of bool array indicating presence/absence of peaks (1) and troughs (2)
flap_sig = np.zeros(len(sig))
flap_sig[peaks] = 1
flap_sig[troughs] = 3

data = np.absolute(sig[troughs].values - sig[peaks].values)
ipt = accFn.interpeaktrough(data)

large_inds = np.absolute(sig[peaks].values - sig[troughs].values) > ipt
flap_peaks, flap_troughs = peaks[large_inds],troughs[large_inds]


# convert to flapping indeces (of acc signal)
flap_inds = np.sort(np.concatenate([flap_peaks,flap_troughs]))
# find gaps between flap signals greater than twice the typical flapping frequency
flap_mask = np.zeros(len(sig))
start,end = accFn.find_gaps(flap_inds,2*fs/flap_freq)
for x,y in zip(start,end):
    flap_mask[x:y] = 1

flap_inds = np.array([i for i, x in enumerate(flap_mask) if x == 1])
starts,ends = accFn.find_gaps(flap_inds, fs*bout_gap)
flap_bouts = np.zeros(len(sig))
for x,y in zip(starts,ends):
    flap_bouts[x:y] = 1

fig,ax=plt.subplots()
ax.plot(test.acc.DT,test.acc.Z)
# ax.scatter(test.acc.DT[flap_peaks],
#             test.acc.DZ[flap_peaks],
#             c="red",
#             )
# ax.scatter(test.acc.DT[flap_troughs],
#             test.acc.DZ[flap_troughs],
#             c="green",
#             )
ax.plot(
    test.acc.DT[flap_mask==1],
    test.acc.Z[flap_mask==1],
    'r+',
)
plt.show()

plt.hist(data,10000)
plt.show()

import numpy as np
from scipy.signal import find_peaks, remez, filtfilt, spectrogram
sig = test.acc.DZ
peaks,_ = find_peaks(sig)
troughs,_ = find_peaks(-sig)
while len(peaks) != len(troughs):
    if peaks[0] > troughs[0]:
        troughs = np.delete(troughs,0)
    if peaks[-1] > troughs[-1]:
        peaks = np.delete(peaks,-1)

data = np.absolute(sig[peaks].values - sig[troughs].values)
ipt = accFn.interpeaktrough(data)

from scipy import stats
kde = stats.gaussian_kde(data)
x = np.linspace(data.min(),data.max(),100)
p = kde(x)

import matplotlib.pyplot as plt
plt.plot(p)
plt.show()

ppeaks,_ = find_peaks(p)
pks,_ = find_peaks(-p[ppeaks[0]:ppeaks[1]])


# # if only using estimated flight periods, recalculate peaks and troughs
# peaks,troughs,_ = accFn.peak_trough(sig)
# retain only sufficiently large z displacements
large_inds = np.absolute(sig[peaks].values - sig[troughs].values) > ipt
flap_peaks, flap_troughs = peaks[large_inds],troughs[large_inds]
# convert to flapping indeces (of acc signal)
flap_inds = np.sort(np.concatenate([flap_peaks,flap_troughs]))
# find gaps between flap signals greater than twice the typical flapping frequency
flap_mask = np.zeros(len(sig))
start,end = find_gaps(flap_inds,2*fs/flap_freq)
for x,y in zip(start,end):
    flap_mask[x:y] = 1

flap_inds = np.array([i for i, x in enumerate(flap_mask) if x == 1])
starts,ends = find_gaps(flap_inds, fs*bout_gap)
flap_bouts = np.zeros(len(sig))
for x,y in zip(starts,ends):
    flap_bouts[x:y] = 1


# create alternative of bool array indicating presence/absence of peaks (1) and troughs (2)
flap_sig = np.zeros(len(sig))
flap_sig[peaks] = 1
flap_sig[troughs] = 3




pk,tr,_ = accFn.peak_trough(test.acc.DZ)


from scipy.signal import find_peaks
pk = find_peaks(test.acc.DZ)[0]
test.acc.DZ[test.acc.index[pk]]

import matplotlib.pyplot as plt

plt.plot(test.acc.DT,test.acc.Z,
         linewidth=.5)

plt.scatter(
    test.acc.DT[test.flap_peaks-1],
    test.acc.Z[test.flap_peaks-1],
    marker='+',
    c='green'
    )
plt.scatter(
    test.acc.DT[test.flap_troughs-1],
    test.acc.Z[test.flap_troughs-1],
    marker='+',
    c='red'
    )
plt.show()

# plt.scatter(test.acc.DT[test.flap_bouts == 1],
#             test.acc.DZ[test.flap_bouts == 1],
#             marker='+', c = 'green')
plt.scatter(test.acc.DT[test.flap_peaks-1],
            test.acc.DZ[test.flap_peaks-1],
            marker='+', c = 'green',
            label='Flapping peak')
plt.scatter(test.acc.DT[test.flap_troughs],
            test.acc.DZ[test.flap_troughs],
            marker='+', c = 'red',
            label="Flapping trough")
plt.legend(loc="upper right")
plt.show()

plt.plot(test.acc.DT,test.acc.DZ,
         linewidth=.5)
for f in test.flInds:
    plt.plot(test.acc.DT[f],
                test.acc.DZ[f],
                c = 'green'
                )
plt.show()



test_out = pd.concat(
    [
        test.acc.DT,
        pd.Series((test.EthBeh == "Rest") & (test.flight == 0)),
        pd.Series(test.flight),
        pd.Series(test.flap_bouts),
        pd.Series(test.gps.lat),
        pd.Series(test.gps.lon),
        pd.Series(test.EthBeh),
        ],
    axis=1
    )
test_out.iloc[:,'Speed'] = []
test_out.rename(columns={
    0:'on_water',
    1:'flapping_extrema',
    2:'flapping_bouts',
    3:'lat',
    4:'lon',
    5:'Ethogram'},
                inplace=True
                )

test_out.to_csv('example_output.csv',date_format="%Y-%m-%d %H/%M/%S.%f")

test_out.columns

test.plot_acc_behaviours('DZ',cols=[
    "#53c0bd","#8b35ea","#5ab8fe","#ffd700"]#,"#ce7e00"
    )

test.rollSum()

b = 'pitch'
getattr(test.acc,b)
test.acc[b]

plot_acc_behaviours(test, acc_sig, cols=None, plot_vid_forage: bool = False):

# check for simplified upsampled behaviours and generate if not present
if (not hasattr(test,'upsampled_beh')) and (test.tag_type == 'dvl'):
    test.upsample_behaviours()
    if (not hasattr(test.upsampled_beh,'beh_simple')) and (test.tag_type == 'dvl'):
        test.test_det_beh_agreement()

cats = np.unique(test.EthBeh)
n_cats = len(cats)
if cols is None:
    # generate distinct colours
    cols = distinctipy.get_colors(n_cats)
if len(cols) != n_cats:
    raise ValueError(
        f"Number of provided colours ({len(cols)}) does not match the number of detected behaviours ({len(cats)})"
    )
# generate behaviour-based line collections
inds, arcs, behavs = test.get_lines_from_string_list(
    string_list=test.EthBeh, line_list=test.acc[acc_sig]
)
# assign relevant colours
arc_colours = []
for x in behavs:
    arc_colours.append(list(compress(cols, cats == x))[0])

import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(6.4, 3.2))
# set axes limits manually because Collections do not take part in autoscaling
ax.set_xlim(0, len(test.EthBeh))
# ax.set_ylim(-6, 6)

line_collection = LineCollection(arcs, linewidths=1, colors=arc_colours)

ax.add_collection(line_collection)
# get order of behaviours and sort for legend
beh_order = [behavs.index(x) for x in cats]
# create list of behaviours by order of appearance
sorted_behaviours = [x for _, x in sorted(zip(beh_order, cats))]
sorted_cols = [x for _, x in sorted(zip(beh_order, cols))]
ax.set_ylim(
    np.min(test.acc[acc_sig]), np.max(test.acc[acc_sig])
)
# generate legend objects
proxies = [test.make_proxy(x, linewidth=1) for x in sorted_cols]
# overlay detected foraging if requested
if plot_vid_forage:
    np.unique(test.upsampled_beh.Behaviour)
    c = np.concatenate(
        (
            np.where(test.upsampled_beh.Behaviour == "s"),
            np.where(test.upsampled_beh.Behaviour == "d"),
        ),
        axis=1,
    )
    c.sort(kind="mergesort")
    ax.plot(c.tolist()[0], [1] * c.shape[1], "k*")
    legend_obj = Line2D([], [], color="k", marker="*", linestyle="None")
    proxies.append(legend_obj)
    ax.plot(
        np.where(test.upsampled_beh.beh_simple == "AT")[0],
        [-1] * sum(test.upsampled_beh.beh_simple == "AT"),
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
ax.title.set_text(f"{test.tagname} estimated behaviours")

plt.show()




flinds, flight = flightestimate(
    test.acc.Z,
    test.rolling_freq_sum,
    test.accfs,
    test.acc.DT,
    numPoints=24*4
)

test.flight_est(numPoints=24*4)

(not(hasattr(test,'flight'))) & (test.tag_type=="dvl")

test.flapping()



flap_times = [x for x,s in zip(test.acc.DT,test.flap_bouts) if s == 1]


sig = test.acc.Z
peaks,troughs,_ = peak_trough(sig)
data = np.absolute(sig[peaks].values - sig[troughs].values)
ipt = interpeaktrough(data)

from scipy.stats import gaussian_kde
kde = gaussian_kde(data)
x = np.linspace(data.min(),data.max(),100)
p = kde(x)

peak_trough(test.acc.Z)

interpeaktrough(data)
# if only using estimated flight periods, recalculate peaks and troughs
peaks,troughs,_ = peak_trough(sig)
# retain only sufficiently large z displacements
large_inds = np.absolute(sig[peaks].values - sig[troughs].values) > ipt
# convert to flapping indeces (of acc signal)
flap_inds = np.sort(np.concatenate([peaks[large_inds],troughs[large_inds]]))
# find gaps between flap signals greater than twice the typical flapping frequency
flap_mask = np.zeros(len(sig))
start,end = find_gaps(flap_inds,2*fs/flap_freq)
for x,y in zip(start,end):
flap_mask[x:y] = 1

flap_inds = np.array([i for i, x in enumerate(flap_mask) if x == 1])
starts,ends = find_gaps(flap_inds, fs*bout_gap)
flap_bouts = np.zeros(len(sig))
for x,y in zip(starts,ends):
flap_bouts[x:y] = 1


test.flight_est(24*4)
test.flInds
toEx = test.calculate_thresholds()
test.beh_detect(toEx)

test.plot_acc_behaviours(test.acc.DZ, plot_vid_forage=False)

(not hasattr(test,'upsampled_beh')) and (test.tag_type == 'dvl')
(not hasattr(test.upsampled_beh,'beh_simple')) and (test.tag_type == 'dvl')

import matplotlib.pyplot as plt
plt.plot(
    test.acc.DT,
    test.acc.Z,
)
plt.plot(
    test.acc.DT[test.flap_peaks],
    test.acc.Z[test.flap_peaks],
    linestyle='none',
    marker='.',
    c='green',
    )
plt.plot(
    test.acc.DT[test.flap_troughs],
    test.acc.Z[test.flap_troughs],
    linestyle='none',
    marker='.',
    c='red',
    )
plt.show()

test.acc.DZ[test.flap_bouts == 1]

readAxy(test.filepath, cols='gps')
import pandas as pd
testdf = pd.read_csv(test.filepath,sep='\t',header=0,usecols=['Timestamp','X','Y','Z'])


dtFormat(testdf.Timestamp[0])

import re
out = re.findall('^(.*)\+',testdf.Timestamp[0])[0]

if re.findall('^(.*)\+',testdf.Timestamp[0]):
    print('yes')
else:
    print('no')

re.sub('[^0-9/ .:-]','',testdf.Timestamp[0])

testdf.Timestamp

test.readin()
test.acc

test

import pandas as pd
test = pd.read_csv(file_loc,sep='\t',header=0,usecols=['Timestamp','location-lat',])
test.columns

readAxy(test.filepath,cols="gps")

import matplotlib.pyplot as plt
fig,(ax0,ax1)=plt.subplots(2,1,sharex=True)
ax0.plot(test.acc.DT,test.acc.pitch)
ax1.plot(test.acc.DT,test.acc.ODmn)
ax2 = ax1.twinx()
ax2.plot(test.gps.DT,test.speed,c='r')
plt.show()


import numpy as np

# define an empty ethogram
test.EthBeh = np.array(["Unknown" for _ in range(len(test.acc))])
test.EthBeh[np.where(test.acc.ODmn < 0.2)] = "Rest"
DiveUp = median([np.mean(test.acc.pitch[x]) for x in test.flInds]) + 2 * median(
    [np.var(test.acc.pitch[x]) for x in test.flInds]
)
DiveDown = median([np.min(test.acc.pitch[x]) for x in test.flInds]) - 30
test.EthBeh[np.where(test.flap_bouts == 1)] = "FL"
# find pitch changes
pitpeaks, pittroughs, _ = accFn.peak_trough(test.acc.pitch)
PitUp = np.array(test.acc.pitch[pitpeaks]) - np.array(
    test.acc.pitch[pittroughs]
)
PitDown = np.array(test.acc.pitch[pittroughs]) - np.array(
    test.acc.pitch[pitpeaks]
)
PitUpL = list(compress(pitpeaks, abs(PitDown) > toEx))
PitDownL = list(compress(pittroughs, PitUp > toEx))
PitLarge = sorted(PitUpL + PitUpL)

# DROP THE TAKE-OFF AS THESE WERE NOT DEFINED FROM VIDEO FOOTAGE
# median_mean_flight_pitch = median([np.mean(test.acc.pitch[x]) for x in test.flInds])
# # examine flight period starts and search for large pitch change
# for fl_start,fl_end in zip(test.flap_start,test.flap_end):
#     idx_range = range(max([fl_start - test.accfs*4,1]),min([fl_start + test.accfs*4,len(test.acc)]))
#     if any(list(filter(lambda a: a in set(PitUpL),set(list(idx_range))))):
#         # find where the nearest increase in pitch is within 4 seconds of flight starting
#         TkoStart = min([PitUpL[((PitUpL > (max([fl_start - test.accfs*4,1]))) & (PitUpL < min([max(fl_start) + test.accfs*4,len(test.acc)])))]])
#         # find pitch peak after this increase in pitch
#         if any(test.acc.pitch[min(pitpeaks[pitpeaks > TkoStart]):] <= median_mean_flight_pitch):
#             TkoEnd = min([fl_end[fl_end > min(pitpeaks[pitpeaks > TkoStart]) + np.where(test.acc.pitch[min(pitpeaks[pitpeaks > TkoStart]):] < median_mean_flight_pitch)[0][0]],TkoStart + test.accfs*5])
#             test.EthBeh[TkoStart:TkoEnd] = ['Takeoff'] * (TkoEnd-TkoStart)

# identify foraging
Pitdif = np.where(np.diff(PitLarge) > (23.3 * test.accfs))[0]
PitOutSd = np.zeros(len(test.acc.pitch), dtype=int).tolist()
# combine all foraging within 23.3 seconds
for b in range(len(Pitdif) - 1):
    if b == 0:
        if PitLarge[Pitdif[b]] - PitLarge[0] > test.accfs * 1.5:
            PitOutSd[PitLarge[0] : PitLarge[Pitdif[0]]] = [1] * (
                PitLarge[Pitdif[0]] - PitLarge[0]
            )
    elif PitLarge[Pitdif[b + 1]] - PitLarge[Pitdif[b] + 1] > test.accfs * 1.5:
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
        [x == "FL" for x in test.EthBeh],
        [x > test.ODmFL for x in test.acc.ODmn],
    )
]
all_pass_indeces = [i for i, x in enumerate(all_pass) if x == 3]
for x in all_pass_indeces:
    test.EthBeh[x] = "Forage"

dives = np.zeros(len(test.acc.pitch), dtype=int).tolist()
ForSt, ForEd = test.find_changes(test.EthBeh, "Forage")
# dives will have a significant drop in pitch to start and are followed
# by an increase later as the bird returns to the surface
for fs, fe in zip(ForSt, ForEd):
    # search for dives prior to foraging indication
    if any(test.acc.pitch[max(1, (fs - test.accfs)) : fe] < DiveDown):
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
                test.acc.pitch[dive_from : (dive_from + (test.accfs * 10))]
                > DiveUp
            ):
                dive_to = pitpeaks[
                    np.argmax(
                        (pitpeaks > dive_from)
                        & (test.acc.pitch[pitpeaks] > DiveUp)
                    )
                ]
                dives[dive_from:dive_to] = [1] * (dive_to - dive_from)
    if any(test.acc.pitch[fs:fe] < DiveDown):
        troughsin = pittroughs[((pittroughs >= fs) & (pittroughs <= fe))]
        diveTroughs = troughsin[test.acc.pitch[troughsin] < DiveDown]
        for dt in diveTroughs:
            if any(test.acc.pitch[dt : (dt + (test.accfs * 10))] > DiveUp):
                if not set(
                    [next((i for i in reversed(pitpeaks) if i < dt))]
                ).issubset(set(PitDownL)):
                    continue
                else:
                    dive_to = pitpeaks[
                        np.argmax(
                            (pitpeaks > dt)
                            & (test.acc.pitch[pitpeaks] > DiveUp)
                        )
                    ]
                    dives[dt:dive_to] = [1] * (dive_to - dt)
for d in np.where(dives == 1)[0]:
    test.EthBeh[d] = "Dive"
