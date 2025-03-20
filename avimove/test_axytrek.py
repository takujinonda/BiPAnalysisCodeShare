from avimove import birdTag
from avimove.forage_detect.src.utils.loadIn import readAxy, dtFormat
import json

json_file = "avimove/axytrek_config.json"

with open(json_file,'r') as f:
    axy_details = json.load(f)


file_loc = 'H:/My Drive/Data/AxyTrek/2017-09/2017-9_S1.csv'
axy_file_config = {
    "filepath": file_loc,
    "tag_type": "axy",
    "tagname": "2017-9",
    "accfs": 25,
    "long_acc_name":["X","Y","Z"],
    "gps_fixes_per_minute": 1,
}
test=birdTag(**axy_file_config)
test.readin()

toEx = test.calculate_thresholds()
test.beh_detect()

test.rollSum()

flinds, flight = flightestimate(
    test.acc.Z,
    test.rolling_freq_sum,
    test.accfs,
    test.acc.DT,
    24*4
)

test.flight_est(numPoints=24*4)

(not(hasattr(test,'flight'))) & (test.tag_type=="dvl")

test.flapping()

flap_times = [x for x,s in zip(test.acc.DT,test.flap_bouts) if s == 1]

from avimove.forage_detect.src.utils.analyseAcc import *
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








import matplotlib.pyplot as plt
plt.plot(
    test.acc.DT,
    test.acc.DZ,
)
plt.plot(
    test.acc.DT[test.flap_bouts == 1],
    test.acc.DZ[test.flap_bouts == 1],
    linestyle=None,
    marker='.',
    c='green',
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