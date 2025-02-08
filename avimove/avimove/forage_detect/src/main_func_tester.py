from forage_detect.main_func import *

dat = readBIP("https://bipsharedata.s3.ap-northeast-1.amazonaws.com/analysis/cal_wind2/Axy/3a66a690-5d53-4bed-a4a7-4d25a1569c4d/3a66a690-5d53-4bed-a4a7-4d25a1569c4d_std.csv",cols='acc')

FkOshima = [39.400,141.998] # capture site
colRemDat = removeNearCS(dat,[39.400,141.998], 5)

# static acceleration
static = dict(zip(["SY","SX","SZ"],[lowEquiFilt(x, 1, 1.5, 25) for x in [colRemDat.Y, colRemDat.X, colRemDat.Z]]))

# dynamic acceleration
dynamic = dict(zip(["DY","DX","DZ"],[x - y for x,y in zip([colRemDat.Y, colRemDat.X, colRemDat.Z],static.values())]))

# add static and dynamic info
colRemDat = colRemDat.join([pd.DataFrame(static).set_index(colRemDat.index),pd.DataFrame(dynamic).set_index(colRemDat.index)])

# split data across dates
dates = np.unique(colRemDat.DT.dt.date)

#%%
# for date in dates:
dayRemDat = colRemDat.loc[colRemDat.DT.dt.date == dates[1],:]

# check if enough data for analysis (> 2 hours)
# if len(dayRemDat <= (25*2*3600)):
#     continue

f,t,Sxx = hammingSpect(dayRemDat.X, 25)
rollSum = rollingSpecSum(Sxx, f, 3, 5)

# def findMaxThenReduce(sig,fs,minGap=5):
rollTest = rollingSpecSum(Sxx, f, 3, 5)

out = []

# now find round(24 * 20)/no.hours occasions of flight 5 mins apart



# while len(out) < 20:
# np.arange(rollTest.idxmax() - ,rollTest.idxmax())
#%%
import matplotlib.pyplot as plt

plt.plot(dayRemDat.DT[(dayRemDat.DT >= (dayRemDat.DT[0] + np.timedelta64(1960,'ms')))& (dayRemDat.DT <= (dayRemDat.DT.iloc[-1] - np.timedelta64(2,'s')))].values,rollSum.values)
plt.show()

#%%

# add ODBA
dayRemDat['ODBA'] = np.sqrt(dayRemDat[['SY','SX','SZ']].pow(2).sum(axis=1))

plt.plot(dayRemDat.DT.values,dayRemDat.ODBA.rolling(60*25, closed = "both", min_periods = 1).mean().values)
plt.show()

#%%
# Plot the signal
import matplotlib.dates as mdates
fig,(ax1,ax2) = plt.subplots(2)
# plt.subplots(211)
ax1.plot(dayRemDat.DT.values,dayRemDat.DX.values)
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H"))

# plt.xlabel('Time')
ax1.set_xlabel('Time')
ax1.set_ylabel('Dorsoventral acceleration')
# plt.ylabel('Dorsoventral acceleration')




# Plot the spectrogram

powerSpectrum, freqenciesFound, time, imageAxis = ax2.specgram(dayRemDat.DX, Fs=25)

# plt.xlabel('Time')

# plt.ylabel('Frequency')
plt.show()
#%%
kde = KernelDensity(kernel='gaussian').fit([dayRemDat.ODBA.values])
plt.plot(kde)

kde.get_params()

from sklearn.neighbors import KernelDensity
b = stats.gaussian_kde(dayRemDat.X)
b(positions)
plt.plot(stats.gaussian_kde(dayRemDat.X))
plt.show()

rollTest.idxmax()

np.floor(25*60*2.5).astype(int)


