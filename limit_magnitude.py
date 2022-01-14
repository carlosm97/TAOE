# -*- coding: utf-8 -*-
"""
Created on Sun Apr 25 14:13:22 2021

@author: Mar Pérez Sar

GOAL : To compute and plot the limit magnitude for an observation night. 
"""


import os
import pandas as pd
import numpy as np
from matplotlib import rcParams,gridspec
#rcParams['text.usetex'] = True
import matplotlib.pylab as plt
import seaborn as sns
from astropy.time import Time
from astropy.io import fits

# ___________________________________DATA_____________________________________
os.chdir('/home/carlos/Desktop/MSc/Segundo/TAOE/Observations/aaxs/tables_topcat') # Poner el path de la carpeta en la que guardéis los datos
# LECTURA DE DATOS
match_estrellas_buenas, match = 'good_stars_5nov_g.csv','stars_5nov_g.csv'
data_gs=pd.read_csv(match_estrellas_buenas,delimiter= ',', index_col=False,header=0)
data_match=pd.read_csv(match,delimiter= ',', index_col=False,header=0)

if match[11]=='_':
    dat = match[6:8]+' '+match[8:11]+' '+match[-5]
else:
     dat = match[6]+' '+match[7:10]+' '+match[-5]   

# CALCULO DE LA MAGNITUD LIMITE
ZP=np.mean(data_gs['dm'])
uZP=np.std(data_gs['dm'])



SNR=data_match['FLUX_AUTO']/data_match['FLUXERR_AUTO']

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
print('limiting magnitude:',Mag_lim)
