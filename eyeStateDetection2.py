# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 14:13:21 2020

@author: oseho
"""
from preprocessing import*


# Alpha
lowI = 8
highI = 13
kI,lI = butter(order,[lowI, highI],"bandpass", False,"ba", fs)
alphaData = lfilter(kI,lI, cleanData,axis=0)  

# Beta
lowII = 13
highII = 30
kII,lII = butter(order,[lowII, highII],"bandpass", False,"ba", fs)
betaData = lfilter(kII,lII, cleanData,axis=0)    

# Delta
lowIII = 0.5
highIII = 4
kIII,lIII = butter(order,[lowIII, highIII],"bandpass", False,"ba", fs)
deltaData = lfilter(kIII,lIII, cleanData,axis=0)   

# Theta
lowIV = 4
highIV = 8
kIV,lIV = butter(order,[lowIV, highIV],"bandpass", False,"ba", fs)
thetaData = lfilter(kIV,lIV, cleanData,axis=0) 

# Gamma
lowV = 30
highV = 58
kV,lV = butter(order,[lowV, highV],"bandpass", False,"ba", fs)
gammaData = lfilter(kV,lV, cleanData,axis=0)  

# PSD deltaData
bandFreq1, deltaPsdFz = signal.welch(deltaData[:,0], fs, nperseg=win)
bandFreq1, deltaPsdCz = signal.welch(deltaData[:,1], fs, nperseg=win)
bandFreq1, deltaPsdPz = signal.welch(deltaData[:,2], fs, nperseg=win)

# PSD thetaData
bandFreq2, thetaPsdFz = signal.welch(thetaData[:,0], fs, nperseg=win)
bandFreq2, thetaPsdCz = signal.welch(thetaData[:,1], fs, nperseg=win)
bandFreq2, thetaPsdPz = signal.welch(thetaData[:,2], fs, nperseg=win)

# PSD betaData
bandFreq3, betaPsdFz = signal.welch(betaData[:,0], fs, nperseg=win) 
bandFreq3, betaPsdCz = signal.welch(betaData[:,1], fs, nperseg=win)
bandFreq3, betaPsdPz = signal.welch(betaData[:,2], fs, nperseg=win)

# PSD alphaData
bandFreq4, alphaPsdFz = signal.welch(alphaData[:,0], fs, nperseg=win)
bandFreq4, alphaPsdCz = signal.welch(alphaData[:,1], fs, nperseg=win)
bandFreq4, alphaPsdPz = signal.welch(alphaData[:,2], fs, nperseg=win)

# PSD gammaData
bandFreq5, gammaPsdFz = signal.welch(gammaData[:,0], fs, nperseg=win) 
bandFreq5, gammaPsdCz = signal.welch(gammaData[:,1], fs, nperseg=win)
bandFreq5, gammaPsdPz = signal.welch(gammaData[:,2], fs, nperseg=win)


dt = 1/sfreq                # dt = 1/sfreq
totalBandPsd = np.array([[alphaPsdCz.size],[alphaPsdFz.size],[alphaPsdPz.size],
                      [betaPsdCz.size],[betaPsdFz.size],[betaPsdPz.size],
                      [deltaPsdCz.size],[deltaPsdFz.size],[deltaPsdPz.size],
                      [thetaPsdCz.size],[thetaPsdFz.size],[thetaPsdPz.size],
                      [gammaPsdCz.size],[gammaPsdFz.size],[gammaPsdPz.size]])
meanBandPsd = np.mean(totalBandPsd)
end = meanBandPsd/sfreq # stop = 155000/sfreq = 310.008s
ts = np.arange(0,end,dt)

# Channel Fz plot
plt.plot(ts,alphaPsdFz, label='alpha') 
plt.plot(ts,gammaPsdFz, label='gamma')
plt.plot(ts,betaPsdFz, label='beta')
plt.plot(ts,deltaPsdFz, label='delta')
plt.plot(ts,thetaPsdFz, label='theta')
plt.xlabel('time in seconds')
plt.ylabel('PSD')
plt.title('Fz: low Frequencies-Gamma PSD Plot')
plt.xlim(0,0.6)
#plt.ylim(0,50)
plt.legend()
plt.show()

# Channel Cz plot
plt.plot(ts,alphaPsdCz, label='alpha') 
plt.plot(ts,gammaPsdCz, label='gamma')
plt.plot(ts,betaPsdCz, label='beta')
plt.plot(ts,deltaPsdCz, label='delta')
plt.plot(ts,thetaPsdCz, label='theta')
plt.xlabel('time in seconds')
plt.ylabel('PSD')
plt.title('Cz: low Frequencies-Gamma PSD Plot')
plt.xlim(0,0.6)
#plt.ylim(0,1000)
plt.legend()
plt.show()

# Channel Pz plot
plt.plot(ts,alphaPsdPz, label='alpha') 
plt.plot(ts,gammaPsdPz, label='gamma')
plt.plot(ts,betaPsdPz, label='beta')
plt.plot(ts,deltaPsdPz, label='delta')
plt.plot(ts,thetaPsdPz, label='theta')
plt.xlabel('time in seconds')
plt.ylabel('PSD')
plt.title('Pz: low Frequencies-Gamma PSD Plot')
plt.xlim(0,0.6)
#plt.ylim(0,1000)
plt.legend()
plt.show()

totalPz = np.concatenate((alphaPsdPz,betaPsdPz,thetaPsdPz,deltaPsdPz,
                          gammaPsdPz), axis=0)
maxPsdPz = np.max(totalPz)
# look for the value of maxPsdPz in each of the bands and extract the index
result = np.where(alphaPsdPz == maxPsdPz) 
idx = result[0]

if idx.size == 0:
    result = np.where(betaPsdPz == maxPsdPz)
idx = result[0]
if idx.size == 0:
    result = np.where(thetaPsdPz == maxPsdPz)
idx = result[0]
if idx.size == 0:
    result = np.where(deltaPsdPz == maxPsdPz)
idx = result[0]
if idx.size == 0:
    result = np.where(gammaPsdPz == maxPsdPz)
idx = result[0]

idx=idx.item() # convert index array to scalar number
domFreq = bandFreq5[idx] # use scalar index number to extract its corresponding value in bandFreq5 array

# check sleep state
if (0.5 <= domFreq <= 4): # range from delta band to beta band
    print("dominant frequency = delta band")
    print("SUBJECT STATE:asleep or not alert")
if (4 <= domFreq <= 8): # range from theta band to beta band
    print("dominant frequency = theta band")
    print("SUBJECT STATE:asleep or not alert")
if (8 <= domFreq <= 13): # range from alpha band to beta band
    print("dominant frequency = alpha band")
    print("SUBJECT STATE:awake or alert")
if (13 <= domFreq <= 30): # range from beta band to beta band
    print("dominant frequency = beta band")
    print("SUBJECT STATE:awake or alert")
if (domFreq > 30):
    print("dominant frequency = gamma band ")
    print("SUBJECT STATE: awake or alert")

































































"""
plt.plot(ts,alphaPsdPz, label='alpha') 
plt.plot(ts,gammaPsdPz, label='gamma')
plt.plot(ts,betaPsdPz, label='beta')
plt.plot(ts,deltaPsdPz, label='delta')
plt.plot(ts,thetaPsdPz, label='theta')
plt.xlabel('time in seconds')
plt.ylabel('PSD')
plt.title('Low Frequencies-Gamma PSD Plot')
plt.xlim(0.16,0.6)
plt.ylim(0,200000)
plt.legend()
plt.show()
"""