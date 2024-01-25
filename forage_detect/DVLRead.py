# DVLRead

from pathlib import Path
import pandas as pd

# parent folder for all DVL data
dvlFolders = "E:/My Drive/PhD/Data/2018Shearwater/DVL"

# acc starttimes
st17008 = pd.to_datetime('31/08/2018 06:00:00', format = "%d/%m/%Y %H:%M:%S")
st18012 = pd.to_datetime('31/08/2018 04:00:00', format = "%d/%m/%Y %H:%M:%S")
st18014 = pd.to_datetime('31/08/2018 06:00:00', format = "%d/%m/%Y %H:%M:%S")
st18017 = pd.to_datetime('31/08/2018 04:00:00', format = "%d/%m/%Y %H:%M:%S")
st18018 = pd.to_datetime('31/08/2018 04:00:00', format = "%d/%m/%Y %H:%M:%S")
starttimes = [st17008,st18012,st18014,st18017,st18018]

# video starttimes
vst17008 = pd.to_datetime('31/08/2018 12:00:00', format = "%d/%m/%Y %H:%M:%S")
vst18012 = pd.to_datetime('31/08/2018 10:00:00', format = "%d/%m/%Y %H:%M:%S")
vst18014 = pd.to_datetime('31/08/2018 12:00:00', format = "%d/%m/%Y %H:%M:%S")
vst18017 = pd.to_datetime('31/08/2018 11:00:00', format = "%d/%m/%Y %H:%M:%S")
vst18018 = pd.to_datetime('31/08/2018 12:00:00', format = "%d/%m/%Y %H:%M:%S")
vidtimes = [vst17008,vst18012,vst18014,vst18017,vst18018]

# times to remove for each section after starttimes
rm17008 = [12.5,12.1,12.4,11.65,11.15]
rm18012 = [12.6,11.8,11.5,10.95,10.35,10.15,9.7,9.5]
rm18014 = [12.65,12.35,11.15,10.8,10.5,10.5,9.6]
rm18017 = [11.3,10.6,10.5,10.15]
rm18018 = [11.65,11.35]
rems = [rm17008,rm18012,rm18014,rm18017,rm18018]

# delay seconds
dl17008 = [x * 60 for x in [0, 21, 32, 57, 87, 120]]
dl18012 = [x * 60 for x in [0, 16, 22, 55, 76, 84,110 + 8/60, 114 + 30/60, 120]]
dl18014 = [x * 60 for x in [0, 6 + 10/60, 64, 75, 82, 98 + 15/60,110, 120]]
dl18017 = [x * 60 for x in [0, 74 + 13/60, 78 + 48/60, 79 + 10/60, 120]]
dl18018 = [x * 60 for x in [0, 33, 120]]
delays = [dl17008,dl18012,dl18014,dl18017,dl18018]

tagnames = ['17008','18012','18014','18017','18018']

# dictionaries storing alignment information and start times
tagInfo = {
    '17008': {
        'accStart' : st17008,
        'vidStart' : vst17008,
        'remSec' : rm17008,
        'delSec' : dl17008
    },
    '18012': {
        'accStart' : st18012,
        'vidStart' : vst18012,
        'remSec' : rm18012,
        'delSec' : dl18012
    },
    '18014': {
        'accStart' : st18014,
        'vidStart' : vst18014,
        'remSec' : rm18014,
        'delSec' : dl18014
    },
    '18017': {
        'accStart' : st18017,
        'vidStart' : vst18017,
        'remSec' : rm18017,
        'delSec' : dl18017
    },
    '18018': {
        'accStart' : st18018,
        'vidStart' : vst18018,
        'remSec' : rm18018,
        'delSec' : dl18018
    }
}

for path in Path(dvlFolders).rglob('*acc*.txt'):
    print(path.name)

