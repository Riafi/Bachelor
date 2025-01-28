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

beam_345 = pd.read_csv('data_band_both_100GHz_dooff.csv')
item_maj_345 = []
item_min_345 = []
item_r_345 = []
bmaj_345 = beam_345[' bmaj_value'].to_numpy()
bmin_345 = beam_345[' bmin_value'].to_numpy()
flux_345 = beam_345[' flux_value'].to_numpy()
flux_345 = 10**3*flux_345
flux_err_345 =beam_345[' flux_error'].to_numpy()
flux_err_345 = 10**3*flux_err_345
    #calculte major and minor axis in pc
bmaj_345 = np.deg2rad(bmaj_345/3600) * dis
bmin_345 = np.deg2rad(bmin_345/3600) * dis
    # calculate geometric mean of bmaj and bmin
r_hl_345 = np.sqrt(bmaj_345*bmin_345)
for i in range(0,len(bmaj_345)):
    item_maj.append(bmaj_345)
    item_min.append(bmin_345)
    item_r.append(r_hl_345)
f = f'{flux} $\pm$ {flux_err}'
dict = {'$S_{{100GHz}}$' : flux , '$S_{{345GHz}}$' : flux_345, '$b_{{maj,100GHz}}$ ': bmaj, '$b_{{min,100GHz}}$ ': bmin, '$R_{{hl,100GHz}}$ ': r_hl }
df = pd.DataFrame(dict)
   
print(df) 
print(df.to_latex(float_format="{:.3f}".format))
