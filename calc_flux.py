import astropy as ap
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.utils.data import get_pkg_data_filename
import matplotlib.colors as clr
import pandas as pd


# read in the values from the data and convert the flux values as well as the flux error values from Jy to mJy
#beam_values reads in the data from the 100GHz sourceband
beam_values = pd.read_csv('data_band_both_100GHz_dooff.csv')
flux = beam_values[' flux_value'].to_numpy()
flux= 10**6*flux
print (flux)
flux_err = beam_values[' flux_error'].to_numpy()
flux_err = 10**6*flux_err
print (type(flux))
#beam_345 reads in the data from the 345GHz sourceband
beam_345 = pd.read_csv('data_band_both_345GHz_dooff.csv')
print (beam_345)
flux345 = beam_345[' flux_value'].to_numpy()
flux345=10**6*flux345
print (flux345)
flux345_err = beam_345[' flux_error'].to_numpy()
flux345_err = flux345_err*10**6
print (type(flux345))

# calculate the dust and free free emission contributions using formulas from paper(!!!)
S_100 = np.linspace (20,2000, 1000)
S_345_25 = 4* S_100*((345/100)**(-0.1))
S_345_50 = 2* S_100*((345/100)**(-0.1))
S_345_100 =  S_100*((345/100)**(-0.1))

S_345 = np.linspace (20,2000, 1000)
S_100_25 = 4* S_345 *((100/345)**3.5)
S_100_50 = 2* S_345 *((100/345)**3.5)
S_100_100 =  S_345 *((100/345)**3.5)


#plotting the flux values from both bands against each other
plt.errorbar(flux, flux345, xerr= flux_err, yerr=flux345_err, color='midnightblue', marker='+',capsize=2,  linestyle='none')
#plotting the emission contributions for 25%,50% and 100%
plt.plot(S_100,S_345_25, color ='dimgray')
plt.plot(S_100,S_345_50, color='darkgrey')
plt.plot(S_100,S_345_100, color = 'lightgray')
plt.plot(S_345,S_100_25, color = 'dimgray')
plt.plot(S_345,S_100_50, color='darkgrey')
plt.plot(S_345,S_100_100, color = 'lightgray')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$S_{100}$ [$\mu$Jy]')
plt.ylabel('$S_{345}$ [$\mu$Jy]')
plt.show()

#calculating the dust contribution at 345 GHz for the YMCs
S_ff = flux*((345/100)**(-0.1))
S_dust = flux345 -S_ff
print (S_dust)
print ('Dust contribution', S_dust/flux345)


#read in the major and minor axis
maj_345 = beam_345[' bmaj_value'].to_numpy()
min_345 = beam_345[' bmin_value'].to_numpy()
#read in the peak intensity of the regions and their error
peak_345 = beam_345[' peak_value'].to_numpy()
peak_345 = 10**3*peak_345 # converting Jy/beam into mJy/beam
peak_e_345 = beam_345[' peak_error'].to_numpy()
peak_e_345 = 10**3*peak_e_345 #converting Jy/beam into mJy/beam
print(peak_345, '\pm', peak_e_345)

#calculate the peak brightness temperature using formula  from website (!!!) 
T_b = 1.222*10**3*(peak_345/(345**2*maj_345*min_345))
print('brightness Temperature of peak intensity at 345GHz is',T_b)

#calculate gas temperature
T_gas = 11.07/(np.log(1+ 11.07/(T_b+0.195)))
print('gas Temperature at 345GHz is', T_gas)