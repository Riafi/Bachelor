import astropy as ap
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.utils.data import get_pkg_data_filename
import matplotlib.colors as clr
import pandas as pd
import scipy


dis = 44*10**6 # Distance to ngc3256






beam = pd.read_csv('data_band_both_100GHz_dooff.csv')
item_maj = []
item_min = []
item_r = []
bmaj = beam[' bmaj_value'].to_numpy()
bmin = beam[' bmin_value'].to_numpy()
flux = beam[' flux_value'].to_numpy()
flux = 10**3*flux
flux_err =beam[' flux_error'].to_numpy()
flux_err = 10**3*flux_err
    #calculte major and minor axis in pc
bmaj = np.deg2rad(bmaj/3600) * dis
bmin = np.deg2rad(bmin/3600) * dis
    # calculate geometric mean of bmaj and bmin
r_hl = np.sqrt(bmaj*bmin)
for i in range(0,len(bmaj)):
    item_maj.append(bmaj)
    item_min.append(bmin)
    item_r.append(r_hl)
f = f'{flux} $\pm$ {flux_err}'
dict = {'S_{100GHz}' : flux ,'b_{maj} ': bmaj, '$b_{min}$ ': bmin, '$R_{hl}$ ': r_hl}
df = pd.DataFrame(dict)
   
print(df) 
print(df.to_latex(float_format="{:.3f}".format))
