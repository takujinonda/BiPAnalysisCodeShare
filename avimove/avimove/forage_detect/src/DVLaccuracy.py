# Script to assess accuracy of behaviour detection method on DVL tag data

import platform
import pandas as pd
import glob

from typing import Union
from birdTag import birdTag

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

accFiles = glob.glob(dvlFolders + '*/*acc*.txt')

tags = [
    {'name':'17008','accSt':'31/08/2018 06:00:00','vidSt':'31/08/2018 12:00:00'},
    {'name':'18012','accSt':'31/08/2018 04:00:00','vidSt':'31/08/2018 10:00:00'},
    {'name':'18014','accSt':'31/08/2018 06:00:00','vidSt':'31/08/2018 12:00:00'},
    {'name':'18017','accSt':'31/08/2018 04:00:00','vidSt':'31/08/2018 11:00:00'},
    {'name':'18018','accSt':'31/08/2018 04:00:00','vidSt':'31/08/2018 12:00:00'}
]

def make_tag_config(
    name: str,
    fname: str,
    tag_type: str,
    acc_fs: int,
    acc_long_name: Union[list,None] = None,
    accSt: str = None,
    vidSt: str = None
    ):
    """
    Create a config file to feed to birdTag class generator.

    Args
    ----
    name
        Tag name
    fname
        Path to data file
    tag_type
        Type of tag (must be recognisable by birdTag)
    acc_fs
        Acceleration data sampling frequnecy (Hz)
    acc_long_name
        list of names for acceleration signals. Must conform to the following
        order - longitudinal (along body axis), dorsoventral (vertical axis),
        lateral
    accSt
        DateTime of acceleration recording start
    vidSt
        DateTime of video recording start

    Returns
    -------
    Dict of details passable to init birdTag class.
    """
    return {
        'filepath': fname,
        'tag_type': tag_type,
        'tagname': fname,
        'accfs': acc_fs,
        'long_acc_name': acc_long_name,
        'gps_fixes_per_minute': gps_fixes_per_minute,
        'accStart': accSt,
        'vidStart': vidSt
    }

dvls = []
for ind,tag in enumerate(tags):
    dvls.append({'name':tag['name'],'data':birdTag(filepath=glob.glob(dvlFolders+tag['name']+'/*acc*.txt')[0],type='dvl',tagname=tag['name'],accfs=20,long_acc_name=['Y','Z','X'],accStart=tag['accSt'],vidStart=tag['vidSt'])})
    dvls[ind]['data'].readin()
    dvls[ind]['data'].readBeh(behavPath=behav_class_loc)
    # calculate flight estimates
    dvls[ind]['data'].flight_est(10,True)
    # calculate flight pitch changes
    dvls[ind]['data'].calculate_thresholds()
    dvls[ind]['data'].flapping()

tags_config = [{'filepath': glob.glob(dvlFolders+x['name']+'/*acc*.txt')[0],\\

}
            filepath: str,
            type: str,
            tagname: str,
            accfs: int,
            long_acc_name: Union[str,None] = None,
            gps_fixes_per_minute: Union[int,None] = None,
            accStart: Union[np.datetime64,None] = None,
            vidStart: Union[np.datetime64,None] = None,
            *args,**kwargs

dvls[0]['data']

##################################################################################
################## CONVERT BEHAVIOURS INTO REGULAR TIMESCALE (1S) ################
##################################################################################

import datetime as dt

def round_seconds(
    obj: dt.datetime
    ) -> dt.datetime:
    """
    Round a datetime to the nearest second.
    """
    if obj.microsecond >= 500_000:
        obj += dt.timedelta(seconds=1)
    return obj.replace(microsecond=0)

test['rounded_time'] = [round_seconds(x) for x in test.Time]

test = dvls[0]['data'].dvl_beh
start_time = test.Time[0]
end_time = test.Time[0] + timedelta(hours=2)
reg_date_range = pd.date_range(start = start_time, end=end_time, freq="1s")
import numpy as np
reg_beh = np.repeat('not set',len(reg_date_range))
from_ind = 0
for idx,x in enumerate(test.rounded_time[:-2]):
    to_ind = np.argmax(reg_date_range > x) - 1
    reg_beh[from_ind:to_ind] = test.Behaviour.iloc[idx]
    from_ind = to_ind
reg_beh[from_ind:] = test.Behaviour.iloc[-1]

regularised_behaviour = pd.DataFrame(
    {
        'time': reg_date_range,
        'beh': reg_beh
    }
)

regularised_behaviour.groupby('beh').count()


test[test.Behaviour == 'AT']



plt.scatter(test.Time,test.index,c=test.Behaviour)
plt.show()

reg_date_range[0:np.argmax(reg_date_range > test.rounded_time[1])]

reg_date_range[0:256]
reg_date_range[256:280]

reg_date_range[257 - 1]



dvls[0]['data'].flapping()






# calculate pitch differences for all tags
medianPitchChanges = []
for x in dvls:
    medianPitchChanges.append(x['data'].calculate_thresholds())
from statistics import median
toEx = median(medianPitchChanges)

dvls[0]['data'].beh_detect(toEx)

test = dvls[0]['data']
ethbeh = np.array(['Unknown' for _ in range(len(test.acc))])

ethbeh[np.where(test.acc.ODmn < .2)] = "Rest"
DiveUp = median([np.mean(test.acc.pitch[x]) for x in test.flInds]) + 2*median([np.var(test.acc.pitch[x]) for x in test.flInds])
DiveDown = median([np.min(test.acc.pitch[x]) for x in test.flInds]) - 30
test.EthBeh[np.where(test.flap_bouts) == 1] = "FL"
# find pitch changes
import utils.analyseAcc as accFn
pitpeaks,pittroughs,_ = accFn.peak_trough(test.acc.pitch)
PitUp = np.array(test.acc.pitch[pitpeaks]) - np.array(test.acc.pitch[pittroughs])
PitDown = np.array(test.acc.pitch[pittroughs]) - np.array(test.acc.pitch[pitpeaks])
from itertools import compress
PitUpL = list(compress(pitpeaks,abs(PitDown) > toEx))
PitDownL = list(compress(pittroughs,PitUp > toEx))
PitLarge = sorted(PitUpL + PitUpL)

Pitdif = np.where(np.diff(PitLarge) > (23.3 * test.accfs))[0]
PitOutSd = np.zeros(len(test.acc.pitch),dtype=int).tolist()
for b in range(len(Pitdif) - 1):
    if b == 0:
        if PitLarge[Pitdif[b]] - PitLarge[0] > test.accfs*1.5:
            PitOutSd[PitLarge[0]:PitLarge[Pitdif[0]]] = [1] * (PitLarge[Pitdif[0]] - PitLarge[0])
    elif PitLarge[Pitdif[b + 1]] - PitLarge[Pitdif[b] + 1] > test.accfs * 1.5:
        PitOutSd[PitLarge[Pitdif[b] + 1]:PitLarge[Pitdif[b+1]]] = [1] * (PitLarge[Pitdif[b+1]] - PitLarge[Pitdif[b] + 1])
        if b == (len(Pitdif) - 2):
            PitOutSd[PitLarge[Pitdif[-1]+1]:PitLarge[-1]] = [1] * (PitLarge[-1] - PitLarge[Pitdif[-1]])

all_pass = [sum(x) for x in zip(PitOutSd, [x == 'FL' for x in test.EthBeh], [x > test.ODmFL for x in test.acc.ODmn])]
pass1 = 
pass2 = 
all_pass_indeces = [i for i, x in enumerate(all_pass) if x == 3]
for x in all_pass_indeces:
    test.EthBeh[x] = "Forage"

dives = np.zeros(len(test.acc.pitch),dtype=int).tolist()
ForSt, ForEd = test.find_changes(test.EthBeh, 'Forage')
# dives will have a significant drop in pitch to start and are followed
# by an increase later as the bird returns to the surface
for fs,fe in zip(ForSt,ForEd):
    # search for dives prior to foraging indication
    if any(test.acc.pitch[max(1,(fs-test.accfs)):fe] < DiveDown):
        dive_from = next((i for i in reversed(pittroughs) if i < fs))
        # check the dive has a large downward pitch to remove possible
        # take-offs
        if not set([next((i for i in reversed(pitpeaks) if i < dive_from))]).issubset(set(PitDownL)):
            continue
        else:
            # search if pitch following potential dive reaches above
            # normal levels
            if any(test.acc.pitch[dive_from:(dive_from + (test.accfs*10))] > DiveUp):
                dive_to = pitpeaks[np.argmax((pitpeaks > dive_from) & (test.acc.pitch[pitpeaks] > DiveUp))]
                dives[dive_from:dive_to] = [1]*(dive_to - dive_from)
    if any(test.acc.pitch[fs:fe] < DiveDown):
        troughsin = pittroughs[((pittroughs >= fs) & (pittroughs <= fe))]
        diveTroughs = troughsin[test.acc.pitch[troughsin] < DiveDown]
        for dt in diveTroughs:
            if any(test.acc.pitch[dt:(dt + (test.accfs * 10))] > DiveUp):
                if not set([next((i for i in reversed(pitpeaks) if i < dt))]).issubset(set(PitDownL)):
                    continue
                else:
                    dive_to = pitpeaks[np.argmax((pitpeaks > dt) & (test.acc.pitch[pitpeaks] > DiveUp))]
                    dives[dt:dive_to] = [1] * (dive_to - dt)
for d in np.where(dives == 1)[0]:
    test.EthBeh[d] = "Dive"



dvls[0]['data'].plot_acc_behaviours('DZ')

test = dvls[0]['data']


import utils.analyseAcc as accFn

import numpy as np
from itertools import compress

ODlow

dvls[0]['data'].plot_acc_behaviours('DZ')

inds, arcs, behavs = dvls[0]['data'].get_changes_in_string_list(dvls[0]['data'].EthBeh,getattr(dvls[0]['data'].acc, 'DZ'))
behs = np.unique(test)
import distinctipy
cols = distinctipy.get_colors(3)
arc_colours = []
for x in behavs:
    arc_colours.append(list(compress(cols,behs == x))[0])

from matplotlib.collections import LineCollection

def flatten(xss):
    return [x for xs in xss for x in xs]
flat_behavs = flatten(behavs)

x = []
n_cats = 3
for b in behs:
    x.append(np.where(test == b)[0][0])
[b for _, b in sorted(zip(x, behs))]
    

fig, ax = plt.subplots(figsize=(6.4, 3.2))
# set axes limits manually because Collections do not take part in autoscaling
# ax.set_xlim(0, len(dvls[0]['data'].EthBeh))
# ax.set_ylim(-6, 6)

line_collection = LineCollection(arcs, colors = arc_colours)

ax.add_collection(line_collection)
# generate legend objects
# ax.legend()
plt.show()

test = dvls[0]['data'].EthBeh
np.unique(test)
a,b,c = dvls[0]['data'].get_changes_in_string_list(test,dvls[0]['data'].acc.DZ)

dvls[0]['data'].EthBeh = np.array(["Unknown" for _ in range(len(dvls[0]['data'].acc))])
dvls[0]['data'].EthBeh[np.where(dvls[0]['data'].acc.ODmn < .2)] = "Rest"

DiveUp = median([np.mean(dvls[0]['data'].acc.pitch[x]) for x in dvls[0]['data'].flInds]) + 2*median([np.var(dvls[0]['data'].acc.pitch[x]) for x in dvls[0]['data'].flInds])
DiveDown = median([np.min(dvls[0]['data'].acc.pitch[x]) for x in dvls[0]['data'].flInds]) - 30

dvls[0]['data'].EthBeh[np.where(dvls[0]['data'].flap_bouts == 1)] = "FL"

sum(dvls[0]['data'].EthBeh == 'FL')


where_flap = [i for i,x in enumerate(dvls[0]['data'].flap_bouts) if x == 1]
for b in where_flap:
    dvls[0]['data'].EthBeh[b] = "FL"
# find pitch changes
pitpeaks,pittroughs,_ = accFn.peak_trough(dvls[0]['data'].acc.pitch)
PitUp = np.array(dvls[0]['data'].acc.pitch[pitpeaks]) - np.array(dvls[0]['data'].acc.pitch[pittroughs])
PitDown = np.array(dvls[0]['data'].acc.pitch[pittroughs]) - np.array(dvls[0]['data'].acc.pitch[pitpeaks])
PitUpL = list(compress(pittroughs,abs(PitDown) > toEx))
PitDownL = list(compress(pitpeaks,PitUp > toEx))
PitLarge = sorted(PitUpL + PitUpL)

median_mean_flight_pitch = median([np.mean(dvls[0]['data'].acc.pitch[x]) for x in dvls[0]['data'].flInds])

# examine flight period starts and search for large pitch change
for fl_start,fl_end in zip(dvls[0]['data'].flap_start,dvls[0]['data'].flap_end):
    # check if any large pitch changes within ± 4 seconds of flight start
    idx_range = range(max([fl_start - dvls[0]['data'].accfs*4,1]),min([fl_start + dvls[0]['data'].accfs*4,len(dvls[0]['data'].acc)]))
    if any(list(filter(lambda a: a in set(PitUpL),set(list(idx_range))))):
        # find where the nearest increase in pitch is within 4 seconds of flight starting
        TkoStart = min(1,[PitUpL[((PitUpL > (max([fl_start - dvls[0]['data'].accfs*4,1]))) & (PitUpL < min([fl_start + dvls[0]['data'].accfs*4,len(dvls[0]['data'].acc)])))]])
        # find pitch peak after this increase in pitch
        if any(dvls[0]['data'].acc.pitch[min(pitpeaks[pitpeaks > TkoStart]):] <= median_mean_flight_pitch):
            TkoEnd = min([fl_end[fl_end > min(pitpeaks[pitpeaks > TkoStart]) + np.where(dvls[0]['data'].acc.pitch[min(pitpeaks[pitpeaks > TkoStart]):] < median_mean_flight_pitch)[0][0]],TkoStart + dvls[0]['data'].accfs*5])
            dvls[0]['data'].EthBeh[TkoStart:TkoEnd] = ['Takeoff'] * (TkoEnd-TkoStart)

# identify foraging
Pitdif = np.where(np.diff(PitLarge) > (23.3 * dvls[0]['data'].accfs))[0]
PitOutSd = np.zeros(len(dvls[0]['data'].acc.pitch),dtype=int).tolist()
# combine all foraging within 23.3 seconds
for b in range(len(Pitdif) - 1):
    if b == 0:
        if PitLarge[Pitdif[b]] - PitLarge[0] > dvls[0]['data'].accfs*1.5:
            PitOutSd[PitLarge[0]:PitLarge[Pitdif[0]]] = [1] * (PitLarge[Pitdif[0]] - PitLarge[0])
    elif PitLarge[Pitdif[b + 1]] - PitLarge[Pitdif[b] + 1] > dvls[0]['data'].accfs * 1.5:
        PitOutSd[PitLarge[Pitdif[b] + 1]:PitLarge[Pitdif[b+1]]] = [1] * (PitLarge[Pitdif[b+1]] - PitLarge[Pitdif[b] + 1])
        if b == (len(Pitdif) - 2):
            PitOutSd[PitLarge[Pitdif[-1]+1]:PitLarge[-1]] = [1] * (PitLarge[-1] - PitLarge[Pitdif[-1]])

all_pass = [sum(x) for x in zip(PitOutSd, [x == 'FL' for x in dvls[0]['data'].EthBeh], [x > dvls[0]['data'].ODmFL for x in dvls[0]['data'].acc.ODmn])]
all_pass_indeces = [i for i, x in enumerate(all_pass) if x == 3]
for x in all_pass_indeces:
    dvls[0]['data'].EthBeh[x] = "Forage"

np.unique(dvls[0]['data'].EthBeh)

dvls[0]['data'].plot_acc_behaviours('DZ')

cats = np.unique(dvls[0]['data'].EthBeh)
n_cats = len(cats)
# generate distinct colours
cols = distinctipy.get_colors(n_cats, pastel_factor=0.7)
# generate behaviour-based line collections
inds, arcs, behavs = dvls[0]['data'].get_changes_in_string_list(dvls[0]['data'].EthBeh,
                                getattr(dvls[0]['data'].acc, 'DZ'))
# extract n colours from list
potential_colours = ["indigo", "blue", "green", "yellow", "orange", "red"]
# assign relevant colours
arc_colours = potential_colours[:n_cats]

fig, ax = plt.subplots(figsize=(6.4, 3.2))
# set axes limits manually because Collections do not take part in autoscaling
ax.set_xlim(0, len(dvls[0]['data'].EthBeh))
ax.set_ylim(-6, 6)

line_collection = LineCollection(arcs,linewidths = 1, color = arc_colours)


ax.add_collection(line_collection)
# generate legend objects
proxies = [dvls[0]['data'].make_proxy(x,linewidth=1) for x in arc_colours]
# get order of behaviours and sort for legend
flat_behavs = dvls[0]['data'].flatten(behavs)
x = []
for b in cats:
    x.append(np.where(dvls[0]['data'].EthBeh == b)[0][0])
ordered_behavs = [b for _, b in sorted(zip(x, cats))]
ax.legend(proxies, ordered_behavs)

plt.show()


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

test = dvls[0]['data'].EthBeh
is_desired = [x == 'Forage' for x in test]
starts = []
ends = []
if test[0] == 'Forage':
    starts.insert(0,0)
[starts.append(i+1) for i,x in enumerate([(sub2-sub1) == 1 for sub1,sub2 in zip(is_desired,is_desired[1:])]) if x]
[ends.append(i+1) for i, x in enumerate([(sub2 - sub1) == -1 for sub1,sub2 in zip(is_desired, is_desired[1:])]) if x]

[(sub2 - sub1) == -11 for sub1,sub2 in zip(is_desired, is_desired[1:])]


dives = np.zeros(len(dvls[0]['data'].acc.pitch),dtype=int).tolist()
ForSt, ForEd = dvls[0]['data'].find_changes(dvls[0]['data'].EthBeh, 'Forage')

list(range(ForSt[0],ForEd[0]))
test[5530:6865]

# dives will have a significant drop in pitch to start and are followed
# by an increase later as the bird returns to the surface
for fs,fe in zip(ForSt,ForEd):
    # search for dives prior to foraging indication
    if any(dvls[0]['data'].acc.pitch[max(1,(fs-dvls[0]['data'].accfs)):fe] < DiveDown):
        dive_from = next((i for i in reversed(pittroughs) if i < fs))
        # check the dive has a large downward pitch to remove possible
        # take-offs
        if not set([next((i for i in reversed(pitpeaks) if i < dive_from))]).issubset(set(PitDownL)):
            continue
        else:
            # search if pitch following potential dive reaches above
            # normal levels
            if any(dvls[0]['data'].acc.pitch[dive_from:(dive_from + (dvls[0]['data'].accfs*10))] > DiveUp):
                dive_to = pitpeaks[np.argmax((pitpeaks > dive_from) & (dvls[0]['data'].acc.pitch[pitpeaks] > DiveUp))]
                dives[dive_from:dive_to] = [1]*(dive_to - dive_from)
    if any(dvls[0]['data'].acc.pitch[fs:fe] < DiveDown):
        troughsin = pittroughs[((pittroughs >= fs) & (pittroughs <= fe))]
        diveTroughs = troughsin[dvls[0]['data'].acc.pitch[troughsin] < DiveDown]
        for dt in diveTroughs:
            if any(dvls[0]['data'].acc.pitch[dt:(dt + (dvls[0]['data'].accfs * 10))] > DiveUp):
                if not set([next((i for i in reversed(pitpeaks) if i < dt))]).issubset(set(PitDownL)):
                    continue
                else:
                    dive_to = pitpeaks[np.argmax((pitpeaks > dt) & (dvls[0]['data'].acc.pitch[pitpeaks] > DiveUp))]
                    dives[dt:dive_to] = [1] * (dive_to - dt)
for d in np.where(dives == 1)[0]:
    dvls[0]['data'].EthBeh[d] = "Dive"

ForSt, ForEd = dvls[0]['data'].find_changes(dvls[0]['data'].EthBeh, 'Forage')
if ForSt != []:
    for fs,fe in zip(ForSt, ForEd):
        if sum(np.diff([PitUpL[x] for x in np.where([(x >= fs) & (x <= fe) for x in PitUpL])[0]]) > (dvls[0]['data'].accfs * 2)) < 2:
            dvls[0]['data'].EthBeh[fs:fe] = ["Failed 1"] * (fe - fs)
        if sum(np.diff([PitDownL[x] for x in np.where([(x >= fs) & (x <= fe) for x in PitDownL])[0]]) > (dvls[0]['data'].accfs * 2)) < 2:
            dvls[0]['data'].EthBeh[fs:fe] = ["Failed 1"] * (fe - fs)

[dvls[0]['data'].EthBeh[x] for x in np.where([x == "Forage" for x in dvls[0]['data'].EthBeh])[0]]

[x.replace('Forage','Take-off bout') for x in dvls[0]['data'].EthBeh]
[x.replace('','Unknown') for x in dvls[0]['data'].EthBeh]




[x >= fs for x in PitUpL]


def get_changes_in_string_list(
        string_list : list,
        line_list : list | None = None
        ):
    """
    Convert list of continuous sections (string or other) and convert to
    structure for LineCollection plotting.
    """
    if line_list is not None:
        assert len(string_list) == len(line_list), "string_list and line_list are not of equal length"
        to_use = line_list
    else:
        to_use = string_list
    lines = []
    strings = []
    inds = []
    idx = 0
    while idx < len(string_list) - 1:
        current_val = string_list[idx]
        if len(np.unique(test[idx:])) == 1:
            break
        next_change = np.where([x != current_val for x in string_list[idx:]])[0][0]
        lines.append(np.array([(x,y) for x,y in zip(list(range(idx,idx+next_change)),
                                      [to_use[b] for b in list(range(idx,idx+next_change))])]))
        strings.append([current_val])
        inds.append(list(range(idx,idx+next_change)))
        idx += next_change

    return inds, lines, strings

from matplotlib.lines import Line2D
def make_proxy(color, **kwargs):
    return Line2D([0, 1], [0, 1], color=color, **kwargs)

import distinctipy
from itertools import compress
def plot_acc_behaviours(
        acc_sig,
        ethogram
        ):

    cats = np.unique(ethogram)
    n_cats = len(cats)
    # generate distinct colours
    cols = distinctipy.get_colors(n_cats)
    behs = np.unique(test)
    # generate behaviour-based line collections
    inds, arcs, behavs = get_changes_in_string_list(ethogram,
                                      acc_sig)
    # assign relevant colours
    arc_colours = []
    for x in behavs:
        arc_colours.append(list(compress(cols,behs == x))[0])

    fig, ax = plt.subplots(figsize=(6.4, 3.2))
    # set axes limits manually because Collections do not take part in autoscaling
    ax.set_xlim(0, len(dvls[0]['data'].EthBeh))
    ax.set_ylim(-6, 6)

    line_collection = LineCollection(arcs,linewidths = 1, colors = arc_colours)

    ax.add_collection(line_collection)
    # generate legend objects
    proxies = [make_proxy(x,linewidth=1) for x in arc_colours]
    ax.legend(proxies, behs)

    plt.show()

plot_acc_behaviours(dvls[0]['data'].acc.DZ, dvls[0]['data'].EthBeh)

test.EthBeh == 'Unknown'


fig, ax = plt.subplots(figsize=(6.4, 3.2))
# set axes limits manually because Collections do not take part in autoscaling
ax.set_xlim(0, len(dvls[0]['data'].EthBeh))
ax.set_ylim(-6, 6)

line_collection = LineCollection(arcs_try[1],linewidths = 1, colors = arc_colors, label=arcs_try[2])
ax.add_collection(line_collection)
plt.legend(handles=)
plt.show()

arcs_try = get_changes_in_string_list(dvls[0]['data'].EthBeh,
                                      dvls[0]['data'].acc.DZ)











        if sum([x >= fs for x in PitUpL] and [x <= fe for x in PitUpL]) < 2:
            dvls[0]['data'].EthBeh[fs:fe] = ["Failed 1"] * (fe - fs)
        if sum(np.diff(PitDownL[(PitDownL >= fs) & (PitDownL <= fe)]) < (dvls[0]['data'].accfs * 2)) < 2:
            dvls[0]['data'].EthBeh[fs:fe] = ["Failed 1"] * (fe - fs)

[x for x in PitUp if str(x) != 'nan']

test = dvls[0]['data']

np.sign(np.array(test.flap_bouts)).diff(1).eq(1) & np.sign(np.array(test.flap_bouts).eq(1))
df[np.sign(df[0]).diff(1).eq(1) & np.sign(df[0]).eq(1)]


transition_points = []
transition_points.append(np.where(test.flap_bouts == 1)[0][0])
# transition back to 0
next = np.where(test.flap_bouts[transition_points[-1]:] == 0)[0][0]
test.flap_bouts[transition_points[-1]:transition_points[-1] + 28]


min([fl_end[fl_end > min(dvls[0].pitpeaks[dvls[0].pitpeaks > TKoStart]) + np.where(dvls[0].pitch[min(dvls[0].pitpeaks[dvls[0].pitpeaks > TKoStart]):] < median_mean_flight_pitch)[0][0]],TKoStart + dvls[0].accfs*5])


test.flap_end[test.flap_end > min(test.pitpeaks)]

np.where(test.flap_bouts[(transition_points[-1]+next):] == 1)[0][0]

def findTransitionPoint(arr): 
    transition_points = []
    # define the first point
    transition_points.append(np.where(arr == 1)[0][0])
    continue_searching = True
    while sum(np.diff(arr[transition_points[-1]:]) != 0):
        next = np.where(arr[transition_points[-1]:] == 0)[0][0]
        transition_points.append(np.where(arr[(transition_points[-1]+next):] == 1)[0][0])
    return transition_points

findTransitionPoint(test.flap_bouts)


    # Initialise lower and upper 
    # bounds 
    lb = 0
    ub = n - 1
  
    # Perform Binary search 
    while (lb <= ub): 
        # Find mid 
        mid = (int)((lb + ub) / 2) 
  
        # update lower_bound if 
        # mid contains 0 
        if (arr[mid] == 0): 
            lb = mid + 1
  
        # If mid contains 1 
        elif (arr[mid] == 1): 
              
            # Check if it is the  
            # left most 1 Return 
            # mid, if yes 
            if (mid == 0 or (mid > 0 and arr[mid - 1] == 0)): 
                return mid 
  
            # Else update  
            # upper_bound 
            ub = mid-1
      
    return -1
  
# Driver code 
arr = [0, 0, 0, 0, 1, 1] 
n = len(arr) 
point = findTransitionPoint(test.flap_bouts, len(test.flap_bouts)); 
if(point >= 0): 
    print("Transition point is ", point) 
else: 
    print("There is no transition point") 
  