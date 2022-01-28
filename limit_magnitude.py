# -*- coding: utf-8 -*-
"""
Created on Sun Apr 25 14:13:22 2021

@author: Mar Pérez Sar

GOAL : To compute and plot the limit magnitude for an observation night. 
"""


import os
import pandas as pd
import numpy as np
#from matplotlib import rcParams,grigspec
#rcParams['text.usetex'] = True
import matplotlib.pylab as plt
import seaborn as sns
from astropy.time import Time
from astropy.io import fits

# ___________________________________DATA_____________________________________
#os.chdir('/home/carlos/Desktop/MSc/Segundo/TAOE/Observations/aaxs/tables_topcat') # Poner el path de la carpeta en la que guardéis los datos

os.chdir('/home/carlos/Desktop/MSc/Segundo/TAOE/Observations/GRUPO1/xfu_table')
#os.chdir('/home/carlos/Desktop/MSc/Segundo/TAOE/Observations/GRUPO1/abqs_tables') 
#os.chdir('/home/carlos/Desktop/MSc/Segundo/TAOE/Observations/aaxs/tables_topcat') 
#os.chdir('/home/carlos/Desktop/MSc/Segundo/TAOE/Observations/zex') 

# LECTURA DE DATOS

match_estrellas_buenas, match = 'SN2021xfu_good_stars_5nov_g.csv','SN2021xfu_stars_5nov_g.csv'
data_gs=pd.read_csv(match_estrellas_buenas,delimiter= ',', index_col=False,header=0)
data_match=pd.read_csv(match,delimiter= ',', index_col=False,header=0)

if match[-12]=='_':
    dat = match[-11:-9]+' '+match[-9:-6]+' '+match[-5]
    dat2 = match[-11:-9]+match[-9:-6]+'_'+match[-5]
else:
     dat = match[-10]+' '+match[-9:-6]+' '+match[-5]   
     dat2 = match[-10]+match[-9:-6]+'_'+match[-5]   

# CALCULO DE LA MAGNITUD LIMITE
ZP=np.mean(data_gs['dm'])
uZP=np.std(data_gs['dm'])
ZP_s = np.mean(data_gs['dm'])-2.5*np.log10(600)
print(ZP,ZP_s,uZP)
#%%
SNR=data_match['FLUX_AUTO']/data_match['FLUXERR_AUTO']
plt.close('all')
cond=((SNR>1) )
plt.scatter(np.log10(SNR[cond]),np.log10(data_match['MAG_APER'][cond]+ZP),color='black',label=dat)
plt.plot((np.linspace(0,4,1000)),np.polyval(np.polyfit(np.log10(SNR[cond]),np.log10(data_match['MAG_APER'][cond]+ZP),1),(np.linspace(0,4,1000))),color='red')
plt.xlabel('log(SNR)',size=17)
plt.ylabel('log(m)',size=17)
plt.legend(fontsize=21,loc='upper right',ncol=1)
plt.grid(which='major', axis='both',alpha=0.3, linestyle='-')
plt.tick_params(length=4, width=0.8, top=False, right=False, labelsize=14)
# plt.xlim(np.min(np.log10(SNR21[cond]))-0.1,np.max(np.log10(SNR21[cond]))+0.1)

# Condicion de magnitud límite: SNR=5.
Mag_lim=10**(np.polyval(np.polyfit(np.log10(data_match['FLUX_AUTO'][cond]/data_match['FLUXERR_AUTO'][cond]),np.log10(data_match['MAG_APER'][cond]+ZP), 1),np.log10(5)))
print('limiting magnitude all stars:',Mag_lim)
SNR=data_gs['FLUX_AUTO']/data_gs['FLUXERR_AUTO']
cond=((SNR>1) )
Mag_lim=10**(np.polyval(np.polyfit(np.log10(data_gs['FLUX_AUTO'][cond]/data_gs['FLUXERR_AUTO'][cond]),np.log10(data_gs['MAG_APER'][cond]+ZP), 1),np.log10(5)))
print('limiting magnitude good stars:',Mag_lim)

plt.savefig('magnitud limite SN zex'+dat+'.png')

#%% -----------------------PLOT OF GOOD STARS----------------------------------
plt.figure(figsize=(8,6))
plt.scatter(data_match[match[-5]+'MeanPSFMag'],data_match['dm'], c='k')
plt.scatter(data_gs[match[-5]+'MeanPSFMag'][data_gs[match[-5]+'MeanPSFMag']>11],data_gs['dm'][data_gs[match[-5]+'MeanPSFMag']>11], c='r',label=dat)
plt.legend(fontsize=21,loc='upper right',ncol=1),plt.grid(which='major', axis='both',alpha=0.3, linestyle='-')
plt.xlabel(match[-5]+'MeanPSFMag',size=17),plt.ylabel('dm',size=17)
plt.xlim(14,19)
plt.ylim(28,29)
plt.tick_params(length=4, width=0.8, top=False, right=False, labelsize=14)
plt.savefig(match_estrellas_buenas[:-3]+'png')


#%%---------------------COMPARACION TIEMPO OSCURO-BRILLANTE--------------------

# NOCHE OSCURA IAC80
os.chdir('/home/carlos/Desktop/MSc/Segundo/TAOE/Observations/GRUPO1/xfu_table')
match_estrellas_buenas, match = 'SN2021xfu_good_stars_2oct_g.csv','SN2021xfu_stars_2oct_g.csv'
data_gs=pd.read_csv(match_estrellas_buenas,delimiter= ',', index_col=False,header=0)
data_match=pd.read_csv(match,delimiter= ',', index_col=False,header=0)

if match[-12]=='_':
    dat = match[-11:-9]+' '+match[-9:-6]+' '+match[-5]
    dat2 = match[-11:-9]+match[-9:-6]+'_'+match[-5]
else:
     dat = match[-10]+' '+match[-9:-6]+' '+match[-5]   
     dat2 = match[-10]+match[-9:-6]+'_'+match[-5]   

# CALCULO DE LA MAGNITUD LIMITE
ZP=np.mean(data_gs['dm'])
uZP=np.std(data_gs['dm'])
ZP_s = np.mean(data_gs['dm'])-2.5*np.log10(600)
print(ZP,ZP_s,uZP)
SNR=data_gs['FLUX_AUTO']/data_gs['FLUXERR_AUTO']
plt.close('all')
cond=((SNR>1) )
plt.scatter(np.log10(SNR[cond]),np.log10(data_gs['MAG_APER'][cond]+ZP),color='black',label='IAC80 dark night observations')
plt.plot((np.linspace(0,4,1000)),np.polyval(np.polyfit(np.log10(SNR[cond]),np.log10(data_gs['MAG_APER'][cond]+ZP),1),(np.linspace(0,4,1000))),color='red',label='IAC80 dark night fit')
plt.xlabel('log(SNR)',size=17)
plt.ylabel('log(m)',size=17)
plt.legend(fontsize=21,loc='upper right',ncol=1)
plt.grid(which='major', axis='both',alpha=0.3, linestyle='-')
plt.tick_params(length=4, width=0.8, top=False, right=False, labelsize=14)
# plt.xlim(np.min(np.log10(SNR21[cond]))-0.1,np.max(np.log10(SNR21[cond]))+0.1)

# Condicion de magnitud límite: SNR=5.
Mag_lim=10**(np.polyval(np.polyfit(np.log10(data_gs['FLUX_AUTO'][cond]/data_gs['FLUXERR_AUTO'][cond]),np.log10(data_gs['MAG_APER'][cond]+ZP), 1),np.log10(5)))
print('limiting magnitude all stars:',Mag_lim)
SNR=data_gs['FLUX_AUTO']/data_gs['FLUXERR_AUTO']
cond=((SNR>1) )
Mag_lim=10**(np.polyval(np.polyfit(np.log10(data_gs['FLUX_AUTO'][cond]/data_gs['FLUXERR_AUTO'][cond]),np.log10(data_gs['MAG_APER'][cond]+ZP), 1),np.log10(5)))
print('limiting magnitude good stars:',Mag_lim)




# NOCHE BRILLANTE INT PARTE OSCURA
os.chdir('/home/carlos/Desktop/MSc/Segundo/TAOE/Observations/GRUPO1/abqs_tables') 
match_estrellas_buenas, match = 'SN2021abqs_good_stars_13nov_g.csv','SN2021abqs_stars_13nov_g.csv'
data_gs=pd.read_csv(match_estrellas_buenas,delimiter= ',', index_col=False,header=0)
data_match=pd.read_csv(match,delimiter= ',', index_col=False,header=0)

if match[-12]=='_':
    dat = match[-11:-9]+' '+match[-9:-6]+' '+match[-5]
    dat2 = match[-11:-9]+match[-9:-6]+'_'+match[-5]
else:
     dat = match[-10]+' '+match[-9:-6]+' '+match[-5]   
     dat2 = match[-10]+match[-9:-6]+'_'+match[-5]   

# CALCULO DE LA MAGNITUD LIMITE
ZP=np.mean(data_gs['dm'])
uZP=np.std(data_gs['dm'])
ZP_s = np.mean(data_gs['dm'])-2.5*np.log10(600)
print(ZP,ZP_s,uZP)

SNR=data_gs['FLUX_AUTO']/data_gs['FLUXERR_AUTO']

cond=((SNR>1) )
plt.scatter(np.log10(SNR[cond]),np.log10(data_gs['MAG_APER'][cond]+ZP),color='blue',label='INT grey night dark observations')
plt.plot((np.linspace(0,4,1000)),np.polyval(np.polyfit(np.log10(SNR[cond]),np.log10(data_gs['MAG_APER'][cond]+ZP),1),(np.linspace(0,4,1000))),color='orange',label='INT grey night dark fit')
plt.xlabel('log(SNR)',size=17)
plt.ylabel('log(m)',size=17)
plt.legend(fontsize=21,loc='upper right',ncol=1)
plt.grid(which='major', axis='both',alpha=0.3, linestyle='-')
plt.tick_params(length=4, width=0.8, top=False, right=False, labelsize=14)
# plt.xlim(np.min(np.log10(SNR21[cond]))-0.1,np.max(np.log10(SNR21[cond]))+0.1)

# Condicion de magnitud límite: SNR=5.
Mag_lim=10**(np.polyval(np.polyfit(np.log10(data_gs['FLUX_AUTO'][cond]/data_gs['FLUXERR_AUTO'][cond]),np.log10(data_gs['MAG_APER'][cond]+ZP), 1),np.log10(5)))
print('limiting magnitude all stars:',Mag_lim)
SNR=data_gs['FLUX_AUTO']/data_gs['FLUXERR_AUTO']
cond=((SNR>1) )
Mag_lim=10**(np.polyval(np.polyfit(np.log10(data_gs['FLUX_AUTO'][cond]/data_gs['FLUXERR_AUTO'][cond]),np.log10(data_gs['MAG_APER'][cond]+ZP), 1),np.log10(5)))
print('limiting magnitude good stars:',Mag_lim)




# NOCHE BRILLANTE INT PARTE BRILLANTE
os.chdir('/home/carlos/Desktop/MSc/Segundo/TAOE/Observations') 
match_estrellas_buenas, match = 'dark_night_gs.txt','dark_night.txt'
data_gs=pd.read_csv(match_estrellas_buenas,delimiter= ',', index_col=False,header=0)
data_match=pd.read_csv(match,delimiter= ',', index_col=False,header=0)

if match[-12]=='_':
    dat = match[-11:-9]+' '+match[-9:-6]+' '+match[-5]
    dat2 = match[-11:-9]+match[-9:-6]+'_'+match[-5]
else:
     dat = match[-10]+' '+match[-9:-6]+' '+match[-5]   
     dat2 = match[-10]+match[-9:-6]+'_'+match[-5]   

# CALCULO DE LA MAGNITUD LIMITE
ZP=np.mean(data_gs['dm'])
uZP=np.std(data_gs['dm'])
ZP_s = np.mean(data_gs['dm'])-2.5*np.log10(600)
print(ZP,ZP_s,uZP)

SNR=data_gs['FLUX_AUTO']/data_gs['FLUXERR_AUTO']

cond=((SNR>1) )
plt.scatter(np.log10(SNR[cond]),np.log10(data_gs['MAG_APER'][cond]+ZP),color='cyan',label='INT grey night dark observations')
plt.plot((np.linspace(0,4,1000)),np.polyval(np.polyfit(np.log10(SNR[cond]),np.log10(data_gs['MAG_APER'][cond]+ZP),1),(np.linspace(0,4,1000))),color='m',label='INT grey night dark fit')
plt.xlabel('log(SNR)',size=17)
plt.ylabel('log(m)',size=17)
plt.legend(fontsize=21,loc='upper right',ncol=1)
plt.grid(which='major', axis='both',alpha=0.3, linestyle='-')
plt.tick_params(length=4, width=0.8, top=False, right=False, labelsize=14)
# plt.xlim(np.min(np.log10(SNR21[cond]))-0.1,np.max(np.log10(SNR21[cond]))+0.1)

# Condicion de magnitud límite: SNR=5.
Mag_lim=10**(np.polyval(np.polyfit(np.log10(data_gs['FLUX_AUTO'][cond]/data_gs['FLUXERR_AUTO'][cond]),np.log10(data_gs['MAG_APER'][cond]+ZP), 1),np.log10(5)))
print('limiting magnitude all stars:',Mag_lim)
SNR=data_gs['FLUX_AUTO']/data_gs['FLUXERR_AUTO']
cond=((SNR>1) )
Mag_lim=10**(np.polyval(np.polyfit(np.log10(data_gs['FLUX_AUTO'][cond]/data_gs['FLUXERR_AUTO'][cond]),np.log10(data_gs['MAG_APER'][cond]+ZP), 1),np.log10(5)))
print('limiting magnitude good stars:',Mag_lim)








plt.legend(fontsize=14)

plt.axvline(x=np.log10(5), ymin=0, ymax=1,linestyle='-.',color='black',alpha=0.6)
plt.text(np.log10(5)-0.1,1.455,'SNR=5', color='black',size=14,rotation=45) 

plt.axvline(x=np.log10(10), ymin=0, ymax=1,linestyle='-.',color='black',alpha=0.6)
plt.text(np.log10(10)-0.1,1.455,'SNR=10', color='black',size=14,rotation=45) 

plt.axvline(x=np.log10(20), ymin=0, ymax=1,linestyle='-.',color='black',alpha=0.6)
plt.text(np.log10(20)-0.1,1.455,'SNR=20', color='black',size=14,rotation=45) 

plt.ylim(1.05,1.45)










