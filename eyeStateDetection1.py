# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 16:43:50 2020

@author: oseho
"""

from defined_modules import*

from eegSignalQuality import*



"""
extract delta band data
filter cleanData to get signals between delta waves (0.5-4Hz)
"""
# band data extraction
lowIII = 0.5
highIII = 4
k1,l1 = butter(order,[lowIII, highIII],"bandpass", False,"ba", fs)
deltaData = lfilter(k1,l1, cleanData,axis=0)         

# conversion of band data to psd format
bandFreq1, deltaPsdFz = signal.welch(deltaData[:,0], fs, nperseg=win)
plt.plot(bandFreq1, deltaPsdFz, label = 'Fz') 
bandFreq1, deltaPsdCz = signal.welch(deltaData[:,1], fs, nperseg=win)
plt.plot(bandFreq1, deltaPsdCz, label = 'Cz')
bandFreq1, deltaPsdPz = signal.welch(deltaData[:,2], fs, nperseg=win)
plt.plot(bandFreq1, deltaPsdPz, label = 'Pz')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Power spectral density (V^2 / Hz)')
plt.title("Delta Band Welch's periodogram")
plt.legend() 
plt.xlim([0, 2.5])
#plt.ylim([0,500])
plt.show()

# marking of band on the psd plot
low, high = 0.5, 4 # Define delta lower and upper limits
# Find intersecting values in frequency vector
idx_deltaFz = np.logical_and(bandFreq1 >= low, bandFreq1 <= high)
idx_deltaCz = np.logical_and(bandFreq1 >= low, bandFreq1 <= high)
idx_deltaPz = np.logical_and(bandFreq1 >= low, bandFreq1 <= high)
# Plot the power spectral density and fill the delta area
plt.plot(bandFreq1, deltaPsdFz)
#plt.plot(bandFreq1, deltaPsdCz)
#plt.plot(bandFreq1, deltaPsdPz)
plt.fill_between(bandFreq1, deltaPsdFz, where=idx_deltaFz, color='skyblue')
#plt.fill_between(bandFreq1, deltaPsdCz, where=idx_deltaCz, color='skyblue')
#plt.fill_between(bandFreq1, deltaPsdPz, where=idx_deltaPz, color='skyblue')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Power spectral density (uV^2 / Hz)')
plt.xlim([0, 6])
plt.title("Delta Band Marking Welch's periodogram")
plt.show()


"""
extract theta band data
filter cleanData to get signals between theta waves (4-8Hz)
"""
# conversion of already extracted thetaData (in eegSignalQuality file)
# to psd format
bandFreq2, thetaPsdFz = signal.welch(thetaData[:,0], fs, nperseg=win)
#plt.plot(bandFreq2, thetaPsdFz, label = 'Fz') 
bandFreq2, thetaPsdCz = signal.welch(thetaData[:,1], fs, nperseg=win)
#plt.plot(bandFreq2, thetaPsdCz, label = 'Cz')
bandFreq2, thetaPsdPz = signal.welch(thetaData[:,2], fs, nperseg=win)
#plt.plot(bandFreq2, thetaPsdPz, label = 'Pz')
#plt.xlabel('Frequency (Hz)')
#plt.ylabel('Power spectral density (V^2 / Hz)')
#plt.title("Theta Band Welch's periodogram")
#plt.legend() 
#plt.xlim([0, 10])
#plt.show()

# marking of band on the psd plot
low, high = 4, 8 # Define delta lower and upper limits
# Find intersecting values in frequency vector
idx_thetaFz = np.logical_and(bandFreq2 >= low, bandFreq2 <= high)
idx_thetaCz = np.logical_and(bandFreq2 >= low, bandFreq2 <= high)
idx_thetaPz = np.logical_and(bandFreq2 >= low, bandFreq2 <= high)
# Plot the power spectral density and fill the delta area
plt.plot(bandFreq2, thetaPsdFz)
plt.plot(bandFreq2, thetaPsdCz)
plt.plot(bandFreq2, thetaPsdPz)
plt.fill_between(bandFreq2, thetaPsdFz, where=idx_thetaFz, color='skyblue')
plt.fill_between(bandFreq2, thetaPsdCz, where=idx_thetaCz, color='red')
plt.fill_between(bandFreq2, thetaPsdPz, where=idx_thetaPz, color='green')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Power spectral density (uV^2 / Hz)')
plt.xlim([0, 10])
plt.title("Theta Band Marking Welch's periodogram")
plt.show()


"""
extract beta band data
filter cleanData to get signals between theta waves (13-30Hz)
"""
# conversion of already extracted betaData (in eegSignalQuality file)
# to psd format
bandFreq3, betaPsdFz = signal.welch(betaData[:,0], fs, nperseg=win)
#plt.plot(bandFreq3, betaPsdFz, label = 'Fz') 
bandFreq3, betaPsdCz = signal.welch(betaData[:,1], fs, nperseg=win)
#plt.plot(bandFreq3, betaPsdCz, label = 'Cz')
bandFreq3, betaPsdPz = signal.welch(betaData[:,2], fs, nperseg=win)
#plt.plot(bandFreq3, betaPsdPz, label = 'Pz')
#plt.xlabel('Frequency (Hz)')
#plt.ylabel('Power spectral density (V^2 / Hz)')
#plt.title("Beta Band Welch's periodogram")
#plt.legend() 
#plt.xlim([0, 33])
#plt.show()

# marking of band on the psd plot
low, high = 13, 30 # Define delta lower and upper limits
# Find intersecting values in frequency vector
idx_betaFz = np.logical_and(bandFreq3 >= low, bandFreq3 <= high)
idx_betaCz = np.logical_and(bandFreq3 >= low, bandFreq3 <= high)
idx_betaPz = np.logical_and(bandFreq3 >= low, bandFreq3 <= high)
# Plot the power spectral density and fill the delta area
plt.plot(bandFreq3, betaPsdFz)
plt.plot(bandFreq3, betaPsdCz)
plt.plot(bandFreq3, betaPsdPz)
plt.fill_between(bandFreq3, betaPsdFz, where=idx_betaFz, color='skyblue')
plt.fill_between(bandFreq3, betaPsdCz, where=idx_betaCz, color='red')
plt.fill_between(bandFreq3, betaPsdPz, where=idx_betaPz, color='green')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Power spectral density (uV^2 / Hz)')
plt.xlim([0, 33])
plt.title("Beta Band Marking Welch's periodogram")
plt.show()


"""
extract alpha band data
filter cleanData to get signals between theta waves (8-13Hz)
"""
# conversion of already extracted alphaData (in eegSignalQuality file)
# to psd format
bandFreq3, alphaPsdFz = signal.welch(alphaData[:,0], fs, nperseg=win)
#plt.plot(bandFreq3, alphaPsdFz, label = 'Fz') 
bandFreq3, alphaPsdCz = signal.welch(alphaData[:,1], fs, nperseg=win)
#plt.plot(bandFreq3, alphaPsdCz, label = 'Cz')
bandFreq3, alphaPsdPz = signal.welch(alphaData[:,2], fs, nperseg=win)
#plt.plot(bandFreq3, alphaPsdPz, label = 'Pz')
#plt.xlabel('Frequency (Hz)')
#plt.ylabel('Power spectral density (V^2 / Hz)')
#plt.title(" Alpha Band Welch's periodogram")
#plt.legend() 
#plt.xlim([0, 15])
#plt.show()

# marking of band on the psd plot
low, high = 8, 13 # Define alpha lower and upper limits
# Find intersecting values in frequency vector
idx_alphaFz = np.logical_and(bandFreq3 >= low, bandFreq3 <= high)
idx_alphaCz = np.logical_and(bandFreq3 >= low, bandFreq3 <= high)
idx_alphaPz = np.logical_and(bandFreq3 >= low, bandFreq3 <= high)
# Plot the power spectral density and fill the delta area
plt.plot(bandFreq3, alphaPsdFz)
plt.plot(bandFreq3, alphaPsdCz)
plt.plot(bandFreq3, alphaPsdPz)
plt.fill_between(bandFreq3, alphaPsdFz, where=idx_alphaFz, color='skyblue')
plt.fill_between(bandFreq3, alphaPsdCz, where=idx_alphaCz, color='red')
plt.fill_between(bandFreq3, alphaPsdPz, where=idx_alphaPz, color='green')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Power spectral density (uV^2 / Hz)')
plt.xlim([0, 15])
plt.title("Alpha Band Marking Welch's periodogram")
plt.show()


"""
extract gamma band data
filter cleanData to get signals between gamma waves (30-58Hz)
"""
# band data extraction
lowIII = 30
highIII = 58
k2,l2 = butter(order,[lowIII, highIII],"bandpass", False,"ba", fs)
gammaData = lfilter(k2,l2, cleanData,axis=0)  
 
# conversion of band data to psd format
bandFreq4, gammaPsdFz = signal.welch(gammaData[:,0], fs, nperseg=win) 
bandFreq4, gammaPsdCz = signal.welch(gammaData[:,1], fs, nperseg=win)
bandFreq4, gammaPsdPz = signal.welch(gammaData[:,2], fs, nperseg=win)


"""
Proof of theories developed in eye state paper 9
"""
dt = 1/sfreq                # dt = 1/sfreq
totalBandPsd = np.array([[alphaPsdCz.size],[alphaPsdFz.size],[alphaPsdPz.size],
                      [betaPsdCz.size],[betaPsdFz.size],[betaPsdPz.size],
                      [deltaPsdCz.size],[deltaPsdFz.size],[deltaPsdPz.size],
                      [thetaPsdCz.size],[thetaPsdFz.size],[thetaPsdPz.size],
                      [gammaPsdCz.size],[gammaPsdFz.size],[gammaPsdPz.size]])
meanBandPsd = np.mean(totalBandPsd)
end = meanBandPsd/sfreq # stop = 155000/sfreq = 310.008s
ts = np.arange(0,end,dt)

# low frequency-High Frequency Power Comparison in the Parietal (Pz) Region
# low frequency = alpha, beta. delta, theta
#  high frequency = gamma and frequencies > 58Hz

# full band plot
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

# alpha-gamma 
plt.plot(ts,alphaPsdPz, label='alpha', color='skyblue') 
plt.plot(ts,gammaPsdPz, label='gamma', color='purple')
plt.xlabel('time in seconds')
plt.ylabel('PSD')
plt.title('Alpha-Gamma PSD Plot')
plt.xlim(0,0.11)
plt.legend()
plt.show()

# beta-gamma
plt.plot(ts,betaPsdPz, label='beta', color='orange')
plt.plot(ts,gammaPsdPz, label='gamma', color='purple')
plt.xlabel('time in seconds')
plt.ylabel('PSD')
plt.title('Beta-Gamma PSD Plot')
plt.xlim(0,0.7)
plt.legend()
plt.show()

# delta-gamma
plt.plot(ts,deltaPsdPz, label='delta', color='green')
plt.plot(ts,gammaPsdPz, label='gamma', color='purple')
plt.xlabel('time in seconds')
plt.ylabel('PSD')
plt.title('Delta-Gamma PSD Plot')
plt.xlim(0,0.07)
plt.legend()
plt.show()

# theta-gamma
plt.plot(ts,thetaPsdPz, label='theta', color='red')
plt.plot(ts,gammaPsdPz, label='gamma', color='purple')
plt.xlabel('time in seconds')
plt.ylabel('PSD')
plt.title('Theta-Gamma PSD Plot')
plt.xlim(0,0.07)
plt.legend()
plt.show()

# low frequency-High Frequency Power Comparison in the Frontal (Fz) Region
        # low frequency = alpha, beta. delta, theta
        #  high frequency = gamma and frequencies > 58Hz
    # alpha-gamma 
plt.plot(ts,alphaPsdFz, label='alpha', color='skyblue') 
plt.plot(ts,gammaPsdFz, label='gamma', color='purple')
plt.xlabel('time in seconds')
plt.ylabel('PSD')
plt.title('Alpha-Gamma PSD Plot')
plt.xlim(0,0.11)
plt.legend()
plt.show()

    # beta-gamma
plt.plot(ts,betaPsdFz, label='beta', color='orange')
plt.plot(ts,gammaPsdFz, label='gamma', color='purple')
plt.xlabel('time in seconds')
plt.ylabel('PSD')
plt.title('Beta-Gamma PSD Plot')
plt.xlim(0,0.7)
plt.legend()
plt.show()

    # delta-gamma
plt.plot(ts,deltaPsdFz, label='delta', color='green')
plt.plot(ts,gammaPsdFz, label='gamma', color='purple')
plt.xlabel('time in seconds')
plt.ylabel('PSD')
plt.title('Delta-Gamma PSD Plot')
plt.xlim(0,0.07)
plt.legend()
plt.show()

    # theta-gamma
plt.plot(ts,thetaPsdFz, label='theta', color='red')
plt.plot(ts,gammaPsdFz, label='gamma', color='purple')
plt.xlabel('time in seconds')
plt.ylabel('PSD')
plt.title('Theta-Gamma PSD Plot')
plt.xlim(0,0.07)
plt.legend()
plt.show()


# Training Data Preparation
# Parietal Region - Pz

    # Eye closed
thetaIdx1 = list(range(1,37))
thetaTrainPz1 = thetaPsdPz[thetaIdx1]
thetaTrainPz1 = thetaTrainPz1.reshape((thetaTrainPz1.size),1)

alphaIdx1 = list(range(1,37)) # index numbers from 30 to 55
alphaTrainPz1 = alphaPsdPz[alphaIdx1] # extract elements with index values from 30 to 55
alphaTrainPz1 = alphaTrainPz1.reshape((alphaTrainPz1.size),1)

betaIdx1 = list(range(1,37))
betaTrainPz1 = betaPsdPz[betaIdx1]
betaTrainPz1 = betaTrainPz1.reshape((betaTrainPz1.size),1)

    # with respect to alpha band
gammaIdx1 = list(range(1,37))  
gammaTrainPz1 = gammaPsdPz[gammaIdx1] #gamma values with respect to alpha band for eye closed
gammaTrainPz1 = gammaTrainPz1.reshape((gammaTrainPz1.size),1)
output10 = ['Eye Closed'] * (gammaTrainPz1.size)
output10 = np.array(output10)
output10 = output10.reshape((output10.size),1)
ecPz1 = np.concatenate((alphaTrainPz1,betaTrainPz1,thetaTrainPz1,
                        gammaTrainPz1,output10), axis=1)

"""
header = ['alpha','beta','theta','gamma','output']
ecPz1Data = pd.DataFrame(ecPz1, columns=header) # convert array into dataframe 
ecPz1Data.to_csv("data1.csv", index=False)# save the dataframe as a csv file 
"""


    # Eye Opened
thetaIdx2 = list(range(130,173))
thetaTrainPz2 = thetaPsdPz[thetaIdx2]
thetaTrainPz2 = thetaTrainPz2.reshape((thetaTrainPz2.size),1)

alphaIdx2 = list(range(130,173)) # index numbers from 30 to 55
alphaTrainPz2 = alphaPsdPz[alphaIdx2] # extract elements with index values from 30 to 55
alphaTrainPz2 = alphaTrainPz2.reshape((alphaTrainPz2.size),1)

betaIdx2 = list(range(130,173))
betaTrainPz2 = betaPsdPz[betaIdx2]
betaTrainPz2 = betaTrainPz2.reshape((betaTrainPz2.size),1)

    # with respect to alpha band
gammaIdx2 = list(range(130,173))  
gammaTrainPz2 = gammaPsdPz[gammaIdx2] #gamma values with respect to alpha band for eye closed
gammaTrainPz2 = gammaTrainPz2.reshape((gammaTrainPz2.size),1)
output11 = ['Eye Opened'] * (gammaTrainPz2.size)
output11 = np.array(output10)
output11 = output10.reshape((output11.size),1)
ecPz2 = np.concatenate((alphaTrainPz2,betaTrainPz2,thetaTrainPz2,
                        gammaTrainPz2,output11), axis=1)
ecPz3 = np.concatenate((ecPz1,ecPz2),axis=0)
 