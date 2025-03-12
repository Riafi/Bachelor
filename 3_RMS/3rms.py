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
beam_values = pd.read_csv('data_band3_3rms_flux_double.csv')
flux = beam_values['flux_value'].to_numpy()
flux= 10**3*flux

flux_err = beam_values['flux_error'].to_numpy()
flux_err = 10**3*flux_err
#beam_345 reads in the data from the 345GHz sourceband
beam_345 = pd.read_csv('data_band7_3rms_flux_double.csv')
flux345 = beam_345['flux_value'].to_numpy()
flux345=10**3*flux345
flux345_err = beam_345['flux_error'].to_numpy()
flux345_err = flux345_err*10**3

#calculating the dust contribution at 345 GHz for the YMCs
S_ff = (flux*10**(-3))*((345/100)**(-0.1))
eS_ff = (flux_err*10**(-3))* ( (345/100)**(-0.1)) #error of free free flux density at 345GHz
S_dust = flux345*10**(-3) -S_ff
eS_dust = np.sqrt((flux345_err*10**(-3))**2 + (eS_ff*10**(-3))**2) # error of flux density at 345GHz using gauss' error estimation
print (' percentage of Dust contribution at 345GHz', (S_dust/flux345)*100)

#calculate free free contribution at 100GHz
S_dust100 = flux345*(100/345)**(3.5)
err_S_dust100 = flux345_err*(100/345)**(3.5)  #error of dust emission flux density at 100GHz
S_ff100 = flux - S_dust100
err_S_ff100 = np.sqrt(flux_err**2 + err_S_dust100**2)      #error of free free emission flux density at 100GHz
print ('dust contribution at 100GHz in flux density', S_dust100)
print ('freefree contribution',S_ff100)
print ('percentage of free free emission at 100GHz', (S_ff100/flux)*100)

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
plt.savefig('S_free_dust_contributions_3rms_double.pdf')
plt.show()

############flux against flux grouped by YMC visibility at diffrent freqencies#########
f_band3_100 =np.take(flux,[0,1,2,3,6,7,8,13,14,15,16])
f_band3_345 =np.take(flux345,[0,1,2,3,6,7,8,13,14,15,16])
f_band3_100_err =np.take(flux_err,[0,1,2,3,6,7,8,13,14,15,16])
f_band3_345_err =np.take(flux345_err,[0,1,2,3,6,7,8,13,14,15,16])
f_band37_100 = np.take(flux,[4,5,12])
f_band37_345 = np.take(flux345,[4,5,12])
f_band37_100_err = np.take(flux_err,[4,5,12])
f_band37_345_err = np.take(flux345_err,[4,5,12])
f_band7_100 =np.take(flux,[9,10,11])
f_band7_345 =np.take(flux345,[9,10,11])
f_band7_100_err =np.take(flux_err,[9,10,11])
f_band7_345_err =np.take(flux345_err,[9,10,11])

print('flux values visible at 3rms at 100',f_band3_100)
fig, ax = plt.subplots()
ax.margins(0)

#plotting the flux values from both bands against each other
plt.errorbar(f_band3_100, f_band3_345, xerr= f_band3_100_err, yerr=f_band3_345_err, color='midnightblue', marker='+',capsize=2,  linestyle='none', zorder = 3.5, label ='100GHz')
plt.errorbar(f_band37_100, f_band37_345, xerr= f_band37_100_err, yerr=f_band37_345_err, color='blueviolet', marker='+',capsize=2,  linestyle='none', zorder = 3.5, label ='100GHz & 345GHz')
plt.errorbar(f_band7_100, f_band7_345, xerr= f_band7_100_err, yerr=f_band7_345_err, color='cornflowerblue', marker='+',capsize=2,  linestyle='none', zorder = 3.5 , label = '345GHz')
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
ax.legend()
plt.savefig('S_free_dust_contributions_3rms_bands_double.pdf')
plt.show()




#calculating the dust contribution at 345 GHz for the YMCs
S_ff = (flux*10**(-3))*((345/100)**(-0.1))
eS_ff = (flux_err*10**(-3))* ( (345/100)**(-0.1)) #error of free free flux density at 345GHz
S_dust = flux345*10**(-3) -S_ff
eS_dust = np.sqrt((flux345_err*10**(-3))**2 + (eS_ff*10**(-3))**2) # error of flux density at 345GHz using gauss' error estimation
print (' percentage of Dust contribution at 345GHz', (S_dust/flux345)*100)

#calculate free free contribution at 100GHz
S_dust100 = flux345*(100/345)**(3.5)
err_S_dust100 = flux345_err*(100/345)**(3.5)  #error of dust emission flux density at 100GHz
S_ff100 = flux - S_dust100
err_S_ff100 = np.sqrt(flux_err**2 + err_S_dust100**2)      #error of free free emission flux density at 100GHz
print ('dust contribution at 100GHz in flux density', S_dust100)
print ('freefree contribution',S_ff100)
print ('percentage of free free emission at 100GHz', (S_ff100/flux)*100)

radi = pd.read_csv('data_band3_3rms _double.csv')
radi3 = pd.read_csv('data_band7_3rms_double.csv')
#read in the major and minor axis
a = radi[' bmaj'].to_numpy()
ar =np.take(a,[0,1,2,3,4,5,6,13,14,15,16])
b = radi3[' bmaj'].to_numpy()
br=b[7:13]
bmaj = np.append(ar,br)
#print (bmaj)
c = radi['bmin'].to_numpy()
cr = np.take(c,[0,1,2,3,4,5,6,13,14,15,16])
d = radi3['bmin'].to_numpy()
dr=d[7:13]
bmin = np.append(cr,dr)
ae = radi[' bmaj_error'].to_numpy()
aer = np.take(ae,[0,1,2,3,4,5,6,13,14,15,16])
be = radi3[' bmaj_error'].to_numpy()
ber=be[7:13]
bmaj_err = np.append(aer,ber)
#print (bmaj)
ce = radi['bmin_error'].to_numpy()
cer = np.take(ce,[0,1,2,3,4,5,6,13,14,15,16])
de = radi3['bmin_error'].to_numpy()
der=de[7:13]
bmin_err = np.append(cer,der)
e = radi['deconvolved_not deconvolved']
f = radi3['deconvolved_not deconvolved']

bmaj_100 = a
bmaj_100_err =ae
bmin_100 = c
bmin_100_err = ce

bmaj_345 = b
bmaj_345_err = be
bmin_345 = d
bmin_345_err = de

print(bmaj,bmin)
#read in the peak intensity of the regions and their error
peak_345 = beam_345['peak_value'].to_numpy()
peak_345 = 10**3*peak_345 # converting Jy/beam into mJy/beam
peak_e_345 = beam_345['peak_error'].to_numpy()
peak_e_345 = 10**3*peak_e_345 #converting Jy/beam into mJy/beam
#print(peak_345, '\pm', peak_e_345)


#calculate gas temperature
T_gas = 11.07/(np.log(1+ 11.07/(45+0.195)))
T_gas_10 = 11.07/(np.log(1+ 11.07/(45+0.195)))*10
T_gas_err = 122.545*0.45/((45+0.195)*(45+11.265)*np.log(11.07/(45+0.195)+1)**2)   #error of gas Temperature
print('gas Temperature at 345GHz is', T_gas)

#calculate dust mass using equation 6 from paper
M_dust= 74220*S_dust*(44**2)*((np.exp(17/T_gas)-1)/0.9)
M_dust_10= 74220*S_dust*(44**2)*((np.exp(17/T_gas_10)-1)/0.9)
M_dust_err =np.sqrt((74220*eS_dust*(44**2)*((np.exp(17/T_gas)-1)/0.9))**2 + (74220/0.9*S_dust*(44**2)*(np.exp(17/T_gas))*17/(T_gas**2)*T_gas_err)**2)   #error of dust mass
print ('dust Mass is', np.log10(M_dust))

#calculating gas mass using equation 7 from paper
M_gas = 120*M_dust
M_gas_10 =120*M_dust_10
M_gas_err = 120*M_dust_err      #error of dust mass
print('gas Mass is', np.log10(M_gas))

#calculate luminosity of the free free emission using S_ff and distance of the ngc3256 at 44 Mpc
L_vt = 4 * np.pi * ((44*10**6*3.1*10**18)**2) * S_ff100*10**(-3)*10**(-23)
L_vt_err = 4*np.pi*(44*10**6*3.1*10**18)**2 * err_S_ff100*10**(-3)*10**(-23)
print ('luminosity ', L_vt) 

#calculate ionizing photon rates using equation 3 from paper sun et al.
Q_0s = 1.6*10**52 * S_ff100 * ((44/10)**2) *((6000/6000)**(-0.51))
Q_0s_err = 1.6*10**52 * err_S_ff100 *(44/10)**2 *((6000/6000)**(-0.51))
print ('Q_0s is', Q_0s)
M_stars = Q_0s / (4.0*10**46)
Mstars_err = Q_0s_err/(4.0*10**46)

print ('stellar mass calculated is', (M_stars))
#calculate ionizing photon rates using equation 9 from paper He et al.
Q_0h = 6.3 * (10**25) * ((1000/1000)**(-0.45))*(100**0.1)*(L_vt)
Q_0h_err = 6.3*(10**25) * ((1000/1000)**(-0.45))*(100*0.1)*(L_vt_err)
M_starh = Q_0h/(4.0*10**46) 
M_starh_err=Q_0h_err/(4.0*10**46)
print ('ionizing photon number and Stellar mass', Q_0h, (M_starh))

print('max and min stellar mass',(min(M_starh)), max(M_starh))
mass=np.linspace(min(M_starh), max(M_starh), 1000)
mas =mass
fig,ax=plt.subplots()
plt.errorbar(M_starh,M_stars,xerr=M_starh_err, yerr=Mstars_err, color='midnightblue', marker='+',capsize=2,  linestyle='none', zorder = 3.5)
plt.plot(mass,mas, color='lightgray',linestyle='--', zorder=2.5)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$M_{stellar}$ [$M_{\odot}$] using He et al. equation')
plt.ylabel('$M_{stellar}$ [$M_{\odot}$] using Sun et al. equation')

plt.savefig('stellar mass comparison_double.pdf')
plt.show()



#total mass
M_tot = M_gas + M_starh
M_tot_err = np.sqrt(M_gas_err**2 + M_starh_err**2)
print(np.log10(M_tot))
gas_fraction = (M_gas/M_tot)*100
gas_fraction_err = np.sqrt((M_gas_err/M_tot)**2 + ((M_gas*M_tot_err)/(M_tot**2))**2)*100



#calculate radii
bmaj_rad = (np.deg2rad(bmaj/3600))/2 * dis
bmaj_rad_err = (np.deg2rad(bmaj_err/3600))/2 * dis
bmin_rad = (np.deg2rad(bmin/3600))/2 * dis
bmin_rad_err = (np.deg2rad(bmin_err/3600))/2 * dis
r_hl = np.sqrt(bmaj_rad*bmin_rad)
r_hl_err = np.sqrt((bmin_rad_err**2*bmaj_rad**2 + bmaj_rad_err**2 * bmin_rad**2)/(bmin_rad*bmaj_rad))
r_hl_plot = r_hl[np.where(r_hl_err/(r_hl*np.log(10))<100)]
r_hl_err_plot = r_hl_err[np.where(r_hl_err/(r_hl*np.log(10))<100)]
M_tot_plot = M_tot[np.where(r_hl_err/(r_hl*np.log(10))<100)]
M_tot_err_plot=M_tot_err[np.where(r_hl_err/(r_hl*np.log(10))<100)]
print ('error of the halflight radius' ,r_hl_err_plot/(r_hl_plot*np.log(10)))
print('not used data in diagramm:',r_hl[np.where(r_hl_err/(r_hl*np.log(10))>100)])
print('beam radii is :', r_hl, r_hl_err)
m=np.linspace(10**4,10**9,1000)
r_beam  = 0*m+ 9.470899419552286
r_beam_better = 0*m +4.73544971
r_beam_half =0*m +0.5*4.73544971/2

r_hl_100 = np.sqrt(((np.deg2rad(bmaj_100/3600))/2 * dis)*((np.deg2rad(bmin_100/3600))/2 * dis))
r_hl_100_err = np.sqrt((((np.deg2rad(bmin_100_err/3600))/2 * dis)**2*((np.deg2rad(bmaj_100/3600))/2 * dis)**2 + ((np.deg2rad(bmaj_100_err/3600))/2 * dis)**2 * ((np.deg2rad(bmin_100/3600))/2 * dis)**2)/(((np.deg2rad(bmin_100/3600))/2 * dis)*((np.deg2rad(bmaj_100/3600))/2 * dis)))

r_hl_345 = np.sqrt(((np.deg2rad(bmaj_345/3600))/2 * dis)*((np.deg2rad(bmin_345/3600))/2 * dis))
r_hl_345_err = np.sqrt((((np.deg2rad(bmin_345_err/3600))/2 * dis)**2*((np.deg2rad(bmaj_345/3600))/2 * dis)**2 + ((np.deg2rad(bmaj_345_err/3600))/2 * dis)**2 * ((np.deg2rad(bmin_345/3600))/2 * dis)**2)/(((np.deg2rad(bmin_345/3600))/2 * dis)*((np.deg2rad(bmaj_345/3600))/2 * dis)))
print('radi of 100GHz', r_hl_100)
print('radi of 345GHz', r_hl_345)

fig,ax = plt.subplots()
ax.margins(0)
plt.errorbar(M_tot_plot,r_hl_plot, yerr = r_hl_err_plot/(r_hl_plot*np.log(10)), xerr= np.abs(M_tot_err_plot/(M_tot_plot*np.log(10))),  color='midnightblue', marker='+',capsize=2,  linestyle='none', zorder = 3.5)
plt.plot(m,r_beam_better, color='lightgray',linestyle='--', zorder=2.5)
plt.text(0.58, 0.6, 'Radius of the restoring beam',color = 'dimgray', rotation_mode = 'default' , rotation = 0, horizontalalignment='left',verticalalignment='top', transform=ax.transAxes)
plt.plot(m,r_beam_half, color='lightgray',linestyle='dotted', zorder=2.5)
plt.text(0.525, 0.225, 'half Radius of the restoring beam',color = 'lightgray', rotation_mode = 'default' , rotation = 0, horizontalalignment='left',verticalalignment='top', transform=ax.transAxes)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$M_{tot}$ [$M_{\odot}$]')
plt.ylabel('$R_{hl}$ [pc]')
ax.set_xlim(left = 5*10**5 , right =2*10**7, auto = True)
ax.set_ylim(bottom = 5*10**(-1) , top = 2*10**1 ,auto =True)
plt.savefig('totalmasstohalflightradius_3rms_double.pdf')
plt.show()
 
################halflight radii to total mass sorted by visiblity in frequency bands#################
r_band3 =np.take(r_hl,[0,1,2,3,6,7,8,13,14,15,16])
M_band3=np.take(M_tot,[0,1,2,3,6,7,8,13,14,15,16])
r_band3err =np.take(r_hl_err,[0,1,2,3,6,7,8,13,14,15,16])
M_band3err =np.take(M_tot_err,[0,1,2,3,6,7,8,13,14,15,16])
r_band37= np.take(r_hl,[4,5,12])
M_band37= np.take(M_tot,[4,5,12])
r_band37err = np.take(r_hl_err,[4,5,12])
M_band37err = np.take(M_tot_err,[4,5,12])
r_band7=np.take(r_hl,[9,10,11])
M_band7=np.take(M_tot,[9,10,11])
r_band7err =np.take(r_hl_err,[9,10,11])
M_band7err =np.take(M_tot_err,[9,10,11])

fig,ax = plt.subplots()
ax.margins(0)
plt.errorbar(M_band3,r_band3, yerr = r_band3err/(r_band3*np.log(10)), xerr= np.abs(M_band3err/(M_band3*np.log(10))),  color='midnightblue', marker='+',capsize=2,  linestyle='none', zorder = 3.5,label ='100GHz')
plt.errorbar(M_band37[np.where(r_band37 > 0.001)],r_band37[np.where(r_band37 > 0.001)], yerr = r_band37err[np.where(r_band37 > 0.001)]/(r_band37[np.where(r_band37 > 0.001)]*np.log(10)), xerr= np.abs(M_band37err[np.where(r_band37 > 0.001)]/(M_band37[np.where(r_band37 > 0.001)]*np.log(10))),  color='blueviolet', marker='+',capsize=2,  linestyle='none', zorder = 3.5,label ='100GHz & 345GHz')
plt.errorbar(M_band7[np.where(r_band7 > 0.001)],r_band7[np.where(r_band7 > 0.001)], yerr = r_band7err[np.where(r_band7 > 0.001)]/(r_band7[np.where(r_band7 > 0.001)]*np.log(10)), xerr= np.abs(M_band7err[np.where(r_band7 > 0.001)]/(M_band7[np.where(r_band7 > 0.001)]*np.log(10))),  color='cornflowerblue', marker='+',capsize=2,  linestyle='none', zorder = 3.5,label ='345GHz')
plt.plot(m,r_beam_better, color='dimgray',linestyle='--', zorder=2.5)
plt.text(0.58, 0.6, 'Radius of the restoring beam',color = 'dimgray', rotation_mode = 'default' , rotation = 0, horizontalalignment='left',verticalalignment='top', transform=ax.transAxes)
plt.plot(m,r_beam_half, color='lightgray',linestyle='dotted', zorder=2.5)
plt.text(0.525, 0.225, 'half Radius of the restoring beam',color = 'lightgray', rotation_mode = 'default' , rotation = 0, horizontalalignment='left',verticalalignment='top', transform=ax.transAxes)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$M_{tot}$ [$M_{\odot}$]')
plt.ylabel('$R_{hl}$ [pc]')
ax.set_xlim(left = 5*10**5 , right =2*10**7, auto = True)
ax.set_ylim(bottom = 5*10**(-1) , top = 2*10**1 ,auto =True)
ax.legend()
plt.savefig('totalmasstohalflightradius_3rms_bands_double.pdf')
plt.show()


ra1 = beam_values[' m0_value'].to_numpy()
D1 = beam_values['m1_value'].to_numpy()
R1 = np.take(ra1,[0,1,2,3,4,5,6,13,14,15,16])
R1 = np.rad2deg(R1)
d1 = np.take(D1,[0,1,2,3,4,5,6,13,14,15,16])
d1 = np.rad2deg(d1)
ra2 = beam_345[' m0_value'].to_numpy()
d2 = beam_345['m1_value'].to_numpy()
R2 = np.rad2deg(ra2[7:13])
d2= np.rad2deg(d2[7:13])
RA = np.append(R1,R2)
dec=np.append(d1,d2)
print(RA,dec)

####gas surface density#########
surface = M_gas/(np.pi*(r_hl**2)*np.log(2))
surface_err = np.sqrt((M_gas_err/(np.pi*(r_hl**2)*np.log(2)))**2 + ((M_gas*2*r_hl_err)/(np.pi*(r_hl**3)*np.log(2)))**2)
print(' gas surface density is', surface)


######optical extinction##########
A_v = 0.0491*surface
A_v_err = 0.0491*surface_err

print('optical extinction of the sources are', A_v)


 
dict = {'$S_{{100GHz}}$' : flux , '$S_{100GHz,ff}$' : S_ff100, '$S_{{345GHz}}$' : flux345, '$S_{345GHz,dust}$' : S_dust , '$M_{Gas}$' : np.log10(M_gas) , '$\Delta M_{gas}$': M_gas_err/(M_gas *np.log(10)),'$R_{hl}$' : r_hl, '$M_{tot}$' : np.log10(M_tot)}
df = pd.DataFrame(dict)
   
print(df) 
print(df.to_latex(float_format="{:.3f}".format, index_names= 'Region ID'))
## flux,flux_err,flux345, flux345_err, S_ff100,err_S_ff100 in mJy
## s_dust,eS_dust in Jy
id = np.linspace(1,17,17)
print(id)
table1 = {'Region_ID':id,'right_Ascension':RA, 'Declination':dec,'flux_100GHz' : flux ,'flux_100GHz_error':flux_err, 'flux_345GHz' : flux345,'flux_345GHz_error': flux345_err,'halflight_Radius':r_hl, 'halflight_Radius_error':r_hl_err }
table1 = pd.DataFrame(table1)
print(table1.to_latex(float_format="{:.3f}".format))

table2 = {'Region_ID':id, 'freefree_flux_100GHz' : S_ff100,'freefree_flux_100GHz_error':err_S_ff100,'dust_flux_345GHz' : S_dust*10**3, 'dust_flux_345GHz_error':eS_dust*10**3, 'gas_Mass':np.log10(M_gas), 'gas_Mass_error':M_gas_err/(M_gas *np.log(10)),'stellar_Mass': np.log10(M_starh),'stellar_Mass_error':M_starh_err/(M_starh *np.log(10)), 'gas_fraction': gas_fraction, 'gas_fraction_error':gas_fraction_err,'total_Mass':np.log10(M_tot), 'total_Mass_error':M_tot_err/(M_tot *np.log(10)), 'halflight_Radius':r_hl, 'halflight_Radius_error':r_hl_err }
table2 = pd.DataFrame(table2)
print(table2.to_latex(float_format="{:.3f}".format))

all_info = {'Region_ID':id, 'right_Ascension':RA, 'Declination':dec,'flux_100GHz' : flux ,'flux_100GHz_error':flux_err, 'freefree_flux_100GHz' : S_ff100,'freefree_flux_100GHz_error':err_S_ff100, 'flux_345GHz' : flux345,'flux_345GHz_error': flux345_err, 'dust_flux_345GHz' : S_dust*10**3, 'dust_flux_345GHz_error':eS_dust*10**3, 'gas_Mass':np.log10(M_gas), 'gas_Mass_error':M_gas_err/(M_gas *np.log(10)),'stellar_Mass': np.log10(M_starh),'stellar_Mass_error':M_starh_err/(M_starh *np.log(10)), 'gas_fraction': gas_fraction, 'gas_fraction_error':gas_fraction_err,'total_Mass':np.log10(M_tot), 'total_Mass_error':M_tot_err/(M_tot *np.log(10)), 'halflight_Radius':r_hl, 'halflight_Radius_error':r_hl_err,'Gas surface density': surface, 'Gas_surface_density_error':surface_err, 'optical_extinction': A_v, 'optical_extinction_error': A_v_err }
for_csv = pd.DataFrame(all_info)
for_csv.to_csv('3_RMS_all_info_double.csv', float_format="{:.3f}".format)

for i in range(0,len(id)):
    print("{:.0f}".format(id[i]), '&', "{:.4f}".format(RA[i]), '&',"{:.4f}".format(dec[i]), '&' , "{:.3f}".format(flux[i]) ,'$\pm$',  "{:.3f}".format(flux_err[i]), '&' ,"{:.3f}".format(flux345[i]) ,'$\pm$',  "{:.3f}".format(flux345_err[i]), '&' ,"{:.1f}".format(r_hl_100[i]) ,'$\pm$' , "{:.1f}".format(r_hl_100_err[i]),'&' ,"{:.1f}".format(r_hl_345[i]) ,'$\pm$' , "{:.1f}".format(r_hl_345_err[i]),'\\\\' )