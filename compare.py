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
radi = pd.read_csv('data_bandboth_100GHz_deconvolved.csv')
radi3 = pd.read_csv('data_bandboth_345GHz_deconvolved.csv')
#read in the major and minor axis
a = radi[' bmaj'].to_numpy()
a =a[0:10]
b = radi3[' bmaj'].to_numpy()
b=b[10:23]
bmaj = np.append(a,b)
print (bmaj)
c = radi['bmin'].to_numpy()
c=c[0:10]
d = radi3['bmin'].to_numpy()
d=d[10:23]
bmin = np.append(c,d)
ae = radi[' bmaj_error'].to_numpy()
ae =ae[0:10]
be = radi3[' bmaj_error'].to_numpy()
be=be[10:23]
bmaj_err = np.append(ae,be)
print (bmaj)
ce = radi['bmin_error'].to_numpy()
ce=ce[0:10]
de = radi3['bmin_error'].to_numpy()
de=de[10:23]
bmin_err = np.append(ce,de)

#calculate radii in pc
bmaj_rad = np.deg2rad(bmaj/3600) * dis
bmaj_rad_err = np.deg2rad(bmaj_err/3600) * dis
bmin_rad = np.deg2rad(bmin/3600) * dis
bmin_rad_err = np.deg2rad(bmin_err/3600) * dis
r_hl = np.sqrt(bmaj_rad*bmin_rad)
r_hl_err = np.sqrt((bmin_rad_err**2*bmaj_rad**2 + bmaj_rad_err**2 * bmin_rad**2)/(bmin_rad*bmaj_rad))
##########################################################
######HWHM for ngc3256##################
dia = 12 #m
f = 345*10**9         #frequency in Hz

beam = np.sqrt(np.deg2rad(0.0452447/3600)* np.deg2rad(0.0435672/3600))
print (beam)
hwhm = beam * dis
print (hwhm)