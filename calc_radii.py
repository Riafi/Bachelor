import astropy as ap
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.utils.data import get_pkg_data_filename
import matplotlib.colors as clr
import pandas as pd
import scipy


dis = 44*10**6 # Distance to ngc3256






beam = pd.read_csv('data_bandboth_deconvolved.csv')
item_maj = []
item_min = []
item_r = []
bmaj = beam[' bmaj'].to_numpy()
bmin = beam[' bmin'].to_numpy()

    #calculte major and minor axis in pc
bmaj = np.deg2rad(bmaj/3600) * dis
bmin = np.deg2rad(bmin/3600) * dis
    # calculate geometric mean of bmaj and bmin
r_hl = np.sqrt(bmaj*bmin)
for i in range(0,len(bmaj)):
    item_maj.append(bmaj)
    item_min.append(bmin)
    item_r.append(r_hl)
decon = beam['converged_deconvolved']
dict = { '$b_{{maj}}$ ': bmaj, '$b_{{min}}$ ': bmin, '$R_{{hl}}$ ': r_hl, 'Region deconvolved' : decon }
df = pd.DataFrame(dict)
   
print(df) 
print(df.to_latex(float_format="{:.3f}".format))
