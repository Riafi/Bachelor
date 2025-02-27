import astropy as ap
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.utils.data import get_pkg_data_filename
import matplotlib.colors as clr
import pandas as pd
import scipy


dis = 44*10**6 # Distance to ngc3256






beam2 = pd.read_csv('data_both_deconvolved_radii.csv')
item_maj = []
item_min = []
item_r = []
bmaj2 = beam2[' bmaj'].to_numpy()
bmaj_err2 = beam2[' bmaj_error'].to_numpy()
bmin2 = beam2['bmin'].to_numpy()
bmin_err2 = beam2['bmin_error'].to_numpy()

    #calculte major and minor axis in pc
bmaj2 = np.deg2rad(bmaj2/3600) * dis
bmaj_err2 = np.deg2rad(bmaj_err2/3600) * dis
bmin2 = np.deg2rad(bmin2/3600) * dis
bmin_err2 = np.deg2rad(bmin_err2/3600) * dis
    # calculate geometric mean of bmaj and bmin
r_hl2 = np.sqrt(bmaj2*bmin2)
r_hl_err2 = np.sqrt((bmin_err2**2*bmaj2**2 + bmaj_err2**2 * bmin2**2)/(bmin2*bmaj2))
for i in range(0,len(bmaj2)):
    item_maj.append(bmaj2)
    item_min.append(bmin2)
    item_r.append(r_hl2)
decon = beam2['deconvolved_not deconvolved']
dict = { '$b_{{maj}}$ ': bmaj2, '$b_{{min}}$ ': bmin2, '$R_{{hl}}$ ': r_hl2, '$\Delta R_{hl}$' : r_hl_err2,  'Region deconvolved' : decon }
df = pd.DataFrame(dict)
   
print(df) 
print(df.to_latex(float_format="{:.3f}".format))

beam = pd.read_csv('data_band3_3rms.csv')
beam3 = pd.read_csv('data_band7_3rms.csv')
bmaj = beam[' bmaj'].to_numpy()
bmaj_err = beam[' bmaj_error'].to_numpy()
bmin = beam['bmin'].to_numpy()
bmin_err = beam['bmin_error'].to_numpy()

bmaj3 = beam3[' bmaj'].to_numpy()
bmaj_err3 = beam3[' bmaj_error'].to_numpy()
bmin3 = beam3['bmin'].to_numpy()
bmin_err3 = beam3['bmin_error'].to_numpy()

    #calculte major and minor axis in pc

bmaj = np.deg2rad(bmaj/3600)/2 * dis
bmaj_err = np.deg2rad(bmaj_err/3600)/2 * dis
bmin = np.deg2rad(bmin/3600)/2 * dis
bmin_err = np.deg2rad(bmin_err/3600)/2 * dis

bmaj3 = np.deg2rad(bmaj3/3600)/2 * dis
bmaj_err3 = np.deg2rad(bmaj_err3/3600)/2 * dis
bmin3 = np.deg2rad(bmin3/3600)/2 * dis
bmin_err3 = np.deg2rad(bmin_err3/3600)/2 * dis

    # calculate geometric mean of bmaj and bmin
r_hl = np.sqrt(bmaj*bmin)
r_hl_err = np.sqrt((bmin_err**2*bmaj**2 + bmaj_err**2 * bmin**2)/(bmin*bmaj))

r_hl3 = np.sqrt(bmaj3*bmin3)
r_hl_err3 = np.sqrt((bmin_err3**2*bmaj3**2 + bmaj_err3**2 * bmin3**2)/(bmin3*bmaj3))

decon = beam['deconvolved_not deconvolved']
decon3=beam3['deconvolved_not deconvolved']
dict1 = { '$R_{{hl,100GHz}}$ ': r_hl, '$\Delta R_{{hl,100GHz}}$ ': r_hl_err,'Region deconvolved at 100GHz' : decon , '$R_{{hl,345GHz}}$ ': r_hl3, '$\Delta R_{hl,345GHz}$' : r_hl_err3, 'Region deconvolved at 345GHz' : decon3 }
df1 = pd.DataFrame(dict1)
   
print(df1) 
print(df1.to_latex(float_format="{:.3f}".format, index_names= 'Region ID'))
df1.to_csv('YMC_Radii.csv', sep='\t', header=True,index=True,index_label='Region ID')
