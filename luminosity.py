import astropy as ap
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.utils.data import get_pkg_data_filename
import matplotlib.colors as clr
import pandas as pd
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score

def gaus(x, a,m,o,d,n):
    return a/(o*(2*np.pi)**(1/2))*np.exp(-(x-m)**2/(2*o)**2) + d*x + n

line = pd.read_csv('/home/riafi/Downloads/member.uid___A001_X2d20_X3a8c.NGC_3256_sci.spw25.cube.regcal.I.pbcor.fits-Z-profile-Region_1-Statistic_Mean-Cooridnate_Current-2025-02-09-05-47-31.tsv', sep='\t',skiprows=5)
velocity = line['# x'].to_numpy()
temp = line['y'].to_numpy()
temp = temp*10**3
temp_rms = temp-9.64541062e-06*velocity +1.03974243e-02

fig,ax= plt.subplots()
plt.plot(velocity,temp,color='darkgray')
plt.plot(velocity,temp,color='midnightblue')
ax.set(ylabel=' flux density per beam [mJy/beam]', xlabel='velocity $V_{lsr}$ /kms$^{-1}$')
plt.show()

fig,ax=plt.subplots()
x0= velocity[np.where((velocity > 2800) & (velocity<2855))]
t0= temp[np.where((velocity> 2800) & (velocity<2855))]
popt, pcov = curve_fit(gaus, x0, t0)
perr = np.sqrt(np.diag(pcov))
print (f'Parameter', popt, perr)
rsquaredfil = r2_score(t0,gaus(x0, *popt))
print (f'$R^2$ ', rsquaredfil)
plt.plot(velocity, temp, label=f'data')
plt.plot(x0, gaus(x0, *popt), color = 'midnightblue', label =f'fit')
ax.set(ylabel=' flux density per beam [mJy/beam]', xlabel='velocity $V_{lsr}$ /kms$^{-1}$')
plt.legend()
plt.show()