import astropy as ap
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.utils.data import get_pkg_data_filename
import matplotlib.colors as clr
import pandas as pd

light = 299792458
k_b=1.380649*10**(-23)
dis=44*10**6

# read in the values from the data and convert the flux values as well as the flux error values from Jy to mJy
#beam_values reads in the data from the 100GHz sourceband
beam_values = pd.read_csv('data_bandboth_100GHz_flux_decon.csv')
flux = beam_values['flux_value'].to_numpy()
flux= 10**3*flux

flux_err = beam_values['flux_error'].to_numpy()
flux_err = 10**3*flux_err
#beam_345 reads in the data from the 345GHz sourceband
beam_345 = pd.read_csv('data_bandboth_345GHz_flux_decon.csv')
flux345 = beam_345['flux_value'].to_numpy()
flux345=10**3*flux345
flux345_err = beam_345['flux_error'].to_numpy()
flux345_err = flux345_err*10**3

# calculate the dust and free free emission contributions using formulas from paper(!!!)
S_100 = np.linspace (0.006,2.5, 1000)

S_345_25 = 4* S_100*((345/100)**(-0.1))
S_345_50 = 2* S_100*((345/100)**(-0.1))
S_345_100 =  S_100*((345/100)**(-0.1))

S_345 = np.linspace (0.1,10, 1000)

S_100_25 = 4* S_345 *((100/345)**3.5)
S_100_50 = 2* S_345 *((100/345)**3.5)
S_100_100 =  S_345*((100/345)**3.5)

fig, ax = plt.subplots()
ax.margins(0)

#plotting the flux values from both bands against each other
plt.errorbar(flux, flux345, xerr= flux_err, yerr=flux345_err, color='midnightblue', marker='+',capsize=2,  linestyle='none', zorder = 3.5)
#plotting the emission contributions for 25%,50% and 100%
plt.plot(S_100,S_345_25, color ='lightgray',zorder=2.5)
plt.text(0.18, 0.1, '$25 \% $',color = 'lightgray', rotation_mode = 'default' , rotation = 45, horizontalalignment='left',verticalalignment='top', transform=ax.transAxes)

plt.plot(S_100,S_345_50, color='darkgrey',zorder=2.5)
plt.text(0.3, 0.1, '$50 \% $', color ='darkgray',  rotation_mode = 'default' , rotation = 45, horizontalalignment='left',verticalalignment='top', transform=ax.transAxes)
plt.plot(S_100,S_345_100, color = 'dimgray',zorder=2.5)
plt.text(0.41, 0.55, '$100 \% $ free-free contribution @ 345GHz',color = 'dimgray', rotation_mode = 'default' , rotation = 44, horizontalalignment='left',verticalalignment='top', transform=ax.transAxes)
plt.plot(S_100_25,S_345,  color = 'lightgray',zorder=2.5)
plt.text(0.025, 0.33, '$25 \% $',color = 'lightgray', rotation_mode = 'default' , rotation = 45, horizontalalignment='left',verticalalignment='top', transform=ax.transAxes)
plt.plot(S_100_50, S_345, color='darkgrey',zorder=2.5)
plt.text(0.025, 0.49, '$50 \% $', color ='darkgray',  rotation_mode = 'default' , rotation = 45, horizontalalignment='left',verticalalignment='top', transform=ax.transAxes)
plt.plot(S_100_100,S_345, color = 'dimgray',zorder=2.5)
plt.text(0, 1, '$100 \% $ dust contribution @ 100GHz',color = 'dimgray', rotation_mode = 'default' , rotation = 44, horizontalalignment='left',verticalalignment='top', transform=ax.transAxes)
ax.fill_between(S_100, S_345_25, S_345_50, color = 'lightgray' , alpha =0.2,zorder = 1.5 )
ax.fill_between(S_100, S_345_50, S_345_100, color = 'dimgray' , alpha =0.2, zorder = 1.5)
ax.fill_between(S_100, S_345_100, 0, color = 'darkgray' , alpha =0.6, zorder = 1.5)
ax.fill_betweenx(S_345, S_100_25, S_100_50, color ='lightgray', alpha = 0.2, zorder = 1.5)
ax.fill_betweenx(S_345, S_100_50, S_100_100, color ='dimgray', alpha = 0.2, zorder = 1.5)
ax.fill_betweenx(S_345, S_100_100, 0, color ='darkgray', alpha = 0.6, zorder = 1.5)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$S_{100}$ [$m$Jy]')
plt.ylabel('$S_{345}$ [$m$Jy]')
ax.set_xlim(left = 0.006 , right =2.5, auto = True)
ax.set_ylim(bottom = 0.05 , top = 5 ,auto =True)
plt.savefig('S_free_dust_contributions.pdf')
plt.show()

#calculating the dust contribution at 345 GHz for the YMCs
S_ff = flux*((345/100)**(-0.1))
eS_ff = flux_err* ( (345/100)**(-0.1)) #error of flux density at 100GHz
S_dust = flux345 -S_ff
eS_dust = np.sqrt(flux345_err**2 + eS_ff**2) # error of flux density at 345GHz using gauss' error estimation
print (' percentage of Dust contribution at 345GHz', (S_dust/flux345)*100)

#calculate free free contribution at 100GHz
S_dust100 = flux345*(100/345)**(3.5)
err_S_dust100 = flux345_err*(100/345)**(3.5)
S_ff100 = flux - S_dust100
err_S_ff100 = np.sqrt(flux_err**2 + err_S_dust100**2)
print ('dust contribution at 100GHz in flux density', S_dust100)
print ('freefree contribution',S_ff100)
print ('percentage of free free emission at 100GHz', (S_ff100/flux)*100)

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
#read in the peak intensity of the regions and their error
peak_345 = beam_345['peak_value'].to_numpy()
peak_345 = 10**3*peak_345 # converting Jy/beam into mJy/beam
peak_e_345 = beam_345['peak_error'].to_numpy()
peak_e_345 = 10**3*peak_e_345 #converting Jy/beam into mJy/beam
#print(peak_345, '\pm', peak_e_345)

#calculate the peak brightness temperature using formula  from website (https://science.nrao.edu/facilities/vla/proposing/TBconv) 
T_b = 1.222*10**3*(peak_345/(345**2*bmaj*bmin))
T_b_err = np.sqrt((1222*peak_e_345/(345**2*bmaj*bmin))**2+((1222*peak_345*bmin_err)/(345**2*bmaj*bmin**2))**2 + ((1222*peak_345*bmaj_err)/(345**2*bmin*bmaj**2))**2)
print('brightness Temperature of peak intensity at 345GHz is',T_b)
#T_B = ((8.69*10**(-4))**2)*S_dust/(2*np.pi*1.381*((np.deg2rad(maj_345/3600)*np.deg2rad(min_345/3600))/4*np.log(2)))
#print ('brightness Temperature,' , T_B)


T= (flux345*10**(-3))*(13.6*(300/345)**2*(1/bmin)*(1/bmaj))
T_err = np.sqrt((flux345_err*10**(-3)*(13.6*(300/345)**2*(1/bmin)*(1/bmaj)))**2 + (flux345*10**(-3)*13.6*(300/345)**2*(bmin_err/(bmin**2))*(1/bmaj))**2 +(flux345*10**(-3)*13.6*(300/345)**2*(bmaj_err/(bmaj**2))*(1/bmin))**2)
print('Temperatur',T)
T_sd = flux345*10**(-3)*0.63*12**2/3514
print ('TemperaturSingledish', T_sd)
#calculate gas temperature
T_gas = 11.07/(np.log(1+ 11.07/(T+0.195)))
T_gas_10 = 11.07/(np.log(1+ 11.07/(T+0.195)))*10
T_gas_err = 122.545*T_err/((T+0.195)*(T+11.265)*np.log(11.07/(T+0.195)+1)**2)
print('gas Temperature at 345GHz is', T_gas)

#calculate dust mass using equation 6 from paper
M_dust= 74.220*S_dust*(44**2)*((np.exp(17/T_gas)-1)/0.9)
M_dust_10= 74.220*S_dust*(44**2)*((np.exp(17/T_gas_10)-1)/0.9)
M_dust_err =np.sqrt((74.220*eS_dust*(44**2)*((np.exp(17/T_gas)-1)/0.9))**2 + (74.220/0.9*S_dust*(44**2)*(np.exp(17/T_gas))*17/(T_gas**2)*T_gas_err)**2)
print ('dust Mass is', np.log10(M_dust))

#calculating gas mass using equation 7 from paper
M_gas = 120*M_dust
M_gas_10 =120*M_dust_10
M_gas_err = 120*M_dust_err
print('gas Mass is', np.log10(M_gas))

#calculate luminosity of the free free emission using S_ff and distance of the ngc3256 at 44 Mpc
L_vt = 4 * np.pi * ((44*10**6)**2) * S_ff100*10**3
print ('luminosity ', L_vt) 

#calculate ionizing photon rates using equation 3 from paper sun et al.
Q_0s = 1.6*10**52 * S_ff100 * ((44/10)**2) *((6000/6000)**(-0.51))
Q_0s_err = 1.6*10**52 * err_S_ff100 *(44/10)**2 *((6000/6000)**(-0.51))
print ('Q_0s is', Q_0s)
M_stars = Q_0s / (4.0*10**46)
Mstars_err = Q_0s_err/(4.0*10**46)

print ('stellar mass calculated is', np.log10(M_stars))
#calculate ionizing photon rates using equation 9 from paper He et al.
Q_0h = 6.3 * (10**25) * ((1000/1000)**(-0.45))*(100**0.1)*(L_vt)
M_starh = Q_0h/(4.0*10**46) 
print ('ionizing photon number and Stellar mass', Q_0h, np.log(M_starh))

#total mass
M_tot = M_gas + M_stars
M_tot_err = np.sqrt(M_gas_err**2 + Mstars_err**2)
print(np.log10(M_tot))
gas_fraction = (M_gas/M_tot)*100


temp_mass = {'$T_b$':T,'$T_{gas}$':T_gas, '$T_{gas,corrected}$':T_gas_10, '$M_{gas}$' :np.log10(M_gas), '$M_{gas,corrected}$':np.log10(M_gas_10)}
df_temp = pd.DataFrame(temp_mass)
print(df_temp.to_latex(float_format="{:.3f}".format, index_names= 'Region ID'))


#calculate radii
bmaj_rad = np.deg2rad(bmaj/3600) * dis
bmaj_rad_err = np.deg2rad(bmaj_err/3600) * dis
bmin_rad = np.deg2rad(bmin/3600) * dis
bmin_rad_err = np.deg2rad(bmin_err/3600) * dis
r_hl = np.sqrt(bmaj_rad*bmin_rad)
r_hl_err = np.sqrt((bmin_rad_err**2*bmaj_rad**2 + bmaj_rad_err**2 * bmin_rad**2)/(bmin_rad*bmaj_rad))
m=np.linspace(10**4,10**9,1000)
r_beam  = 0*m+ 19.198621772
fig,ax = plt.subplots()
ax.margins(0)
plt.errorbar(M_tot,r_hl, yerr = r_hl_err/(np.log(10)), xerr= np.abs(M_tot_err/(np.log(10))),  color='midnightblue', marker='+',capsize=2,  linestyle='none', zorder = 3.5)
plt.plot(m,r_beam, color='lightgray',linestyle='--', zorder=2.5)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$M_{tot}$ [$M_{\odot}$]')
plt.ylabel('$R_{hl}$ [pc]')
ax.set_xlim(left = 3*10**6 , right =4*10**8, auto = True)
ax.set_ylim(bottom = 2*10**(0) , top = 4*10**1 ,auto =True)
plt.savefig('totalmasstohalflightradius.pdf')
plt.show()
 
dict = {'$S_{{100GHz}}$' : flux , '$S_{100GHz,ff}$' : S_ff100, '$S_{{345GHz}}$' : flux345, '$S_{345GHz,dust}$' : S_dust , '$M_{Gas}$' : np.log10(M_gas) , '$\Delta M_{gas}$': M_gas_err/(M_gas *np.log(10)),'$f_{gas}$' : gas_fraction, '$M_{stellar}$' : np.log10(M_stars)}
df = pd.DataFrame(dict)
   
print(df) 
print(df.to_latex(float_format="{:.3f}".format, index_names= 'Region ID'))