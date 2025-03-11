import astropy as ap
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.utils.data import get_pkg_data_filename
import matplotlib.colors as clr
import pandas as pd

#important constants
dis=44*10**6        #distance of ngc3256 from ALMA
light = 299792458          #speed of light

#####radius of my sources#########
radi = pd.read_csv('3_RMS/data_band3_3rms _double.csv')
radi3 = pd.read_csv('3_RMS/data_band7_3rms_double.csv')
#read in the major and minor axis
a = radi[' bmaj'].to_numpy()
ar =np.take(a,[0,1,2,3,4,5,6,7,8,13,14,15,16])
b = radi3[' bmaj'].to_numpy()
b=b[9:13]
bmaj = np.append(ar,b)
#print (bmaj)
c = radi['bmin'].to_numpy()
cr = np.take(c,[0,1,2,3,4,5,6,7,8,13,14,15,16])
d = radi3['bmin'].to_numpy()
d=d[9:13]
bmin = np.append(cr,d)
ae = radi[' bmaj_error'].to_numpy()
aer = np.take(ae,[0,1,2,3,4,5,6,7,8,13,14,15,16])
be = radi3[' bmaj_error'].to_numpy()
be=be[9:13]
bmaj_err = np.append(aer,be)
#print (bmaj)
ce = radi['bmin_error'].to_numpy()
cer = np.take(ce,[0,1,2,3,4,5,6,7,8,13,14,15,16])
de = radi3['bmin_error'].to_numpy()
de=de[9:13]
bmin_err = np.append(cer,de)
#calculate radii
bmaj_rad = (np.deg2rad(bmaj/3600))/2 * dis
bmaj_rad_err = (np.deg2rad(bmaj_err/3600))/2 * dis
bmin_rad = (np.deg2rad(bmin/3600))/2 * dis
bmin_rad_err = (np.deg2rad(bmin_err/3600))/2 * dis
r_hl = np.sqrt(bmaj_rad*bmin_rad)
r_hl_err = np.sqrt((bmin_rad_err**2*bmaj_rad**2 + bmaj_rad_err**2 * bmin_rad**2)/(bmin_rad*bmaj_rad))


print('not used data in diagramm:',r_hl[np.where(r_hl_err/(r_hl*np.log(10))>100)])
print('beam radii is :', r_hl, r_hl_err)

##########################################################
######HWHM for ngc3256##################
dia = 12 #m
f = 345*10**9         #frequency in Hz

fwhm= ((1.02*(light/(f*dia))))*dis/2
print(fwhm)
beam = np.sqrt(np.deg2rad(0.0452447/3600)* np.deg2rad(0.0435672/3600))
print (beam)
hwhm = beam * dis/2
print (hwhm)

###########################################################
########### data sets for comparison ######################
hwhm = 0*r_hl +hwhm                  #HWHM for NGC3256 at R_hl

#sun et al
r_hl_sun = np.array([4.3,2.9,3.7,3.3,1.4,3.5,6.0,4.6,4.1,2.6,4.0,4.5,4.5,2.9,4.0,3.4,4.2,5.9])
r_hl_sun_err = np.array([0.6,1.5,0.4,0.3,3.4,0.6,1.0,0.7,1.2,1.1,1.6,0.7,0.7,1.3,0.6,0.4,1.7,0.9])

hwhm_sun = 0*r_hl_sun+ 8.2/2

#he et al
r_hl_he =np.array([4.5,6.4,9.1,15,7.5,3.2])
r_hl_he_err =np.array([0.5,1.1,2.1,3,1.6,1.1])

hwhm_he = 0*r_hl_he +12/2

#emig et al.
r_hl_emig =np.array([3.0,2.6,3.1,2.5,2.9,2.7,2.4,2.4,2.9,1.4,2.7,2.5,2.2,2.5,2.3,3.1,2.4,3.7,3.1,2.4,2.6,3.1,3.9,3.4,2.2,4.0,2.9])
hwhm_emig = 0*r_hl_emig + 2.2

#levy et al.
r_hl_levy = np.array([0.59,0.19,0.47,0.50,0.76,0.59,1.20,0.46,1.32,0.13,1.29,0.36,0.53])
hwhm_levy = 0*r_hl_levy + 0.48/2
###############################################################
################# diagram plot ################################
hw = np.linspace(0.1, 12, 1000 )
rhl = hw
rhl3 = 3*hw

fig, ax= plt.subplots()
ax.margins(0)
plt.scatter(hwhm, r_hl,  color = 'midnightblue', marker = '+', label = 'NGC3256 (this work)', zorder = 3)
plt.scatter(hwhm_sun, r_hl_sun, color ='cornflowerblue', marker ='o', label = 'NGC3351 (Sun et al. 2024)', zorder = 3)
plt.scatter(hwhm_he, r_hl_he, color = 'blueviolet', marker = 's', label ='Antennae (He et al. 2022)', zorder = 3)
plt.scatter(hwhm_emig, r_hl_emig, color = 'plum', marker ='*', label = 'NGC4945 (Emig et al. 2020)', zorder = 3)
plt.scatter(hwhm_levy, r_hl_levy, color='teal', marker='D', label = 'NGC253 (Levy et al. 2021)', zorder = 3)
plt.plot(hw, rhl, color ='lightgray', linestyle = '--', zorder = 2)
plt.plot(hw, rhl3, color='lightgray', linestyle = 'dotted', zorder = 2)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Beam HWHM [pc]')
plt.ylabel('YMC $R_{hl}$ [pc]')
ax.set_xlim(left = 0.1, right =12, auto = True)
ax.set_ylim(bottom = 0.1 , top = 12 ,auto =True)
ax.legend()
plt.savefig('3_RMS/Compare_HWHM_Radi.pdf')
plt.show()