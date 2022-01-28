#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 24 10:34:06 2022

@author: carlos
"""

# Cargamos el espectro de TNS + templates de Nugent

import os
import pandas as pd
import numpy as np
from matplotlib import gridspec#,rcParams
import matplotlib.pylab as plt
#plt.rcParams['text.usetex'] = True
#import seaborn as sns
from astropy.time import Time
from scipy.interpolate import interp1d
from astropy.io import ascii

os.chdir("/home/carlos/Desktop/MSc/Segundo/TAOE/Observations/SN_LC")
SN = 'aaxs'
if SN == 'aaxs':
    spec = ascii.read('tns_2021aaxs.ascii')
    z = 0.026
    spec['lambda'] /= (1+z) 
elif SN == 'zex':
    spec = ascii.read('tns_2021zex.ascii')
    z = 0.031
    spec['lambda'] /= (1+z)
sn2l = pd.read_csv('sn2l_flux.v1.2.dat',delimiter= '\s+',index_col=False,header=0)
sn2p = pd.read_csv('sn2p_flux.v1.2.dat',delimiter= '\s+',index_col=False,header=0)
sn2n = pd.read_csv('sn2n_flux.v2.1.dat',delimiter= '\s+',index_col=False,header=0)


plt.close('all')
plt.figure(),plt.grid()

spec['flux'] /=np.nanmax(spec['flux'])


plt.plot(spec['lambda'],spec['flux'])
plt.xlim(np.nanmin(spec['lambda']),np.nanmax(spec['lambda']))
0,1,6,11,16,21,36
d = 6

plt.plot(sn2p['lambda'][sn2p['day']==d],\
         sn2p['flux'][sn2p['day']==d]/np.nanmax(sn2p['flux'][sn2p['day']==d])*np.nanmax(spec['flux']),\
         label='SN IIp '+str(d)+' days')


plt.plot(sn2l['lambda'][sn2l['day']==d],\
         sn2l['flux'][sn2l['day']==d]/np.nanmax(sn2l['flux'][sn2l['day']==d])*np.nanmax(spec['flux']),\
         label='SN IIl '+str(d)+' days',alpha=0.5)
plt.ylabel('Relative flux',size=17),plt.xlabel('$\lambda [\AA]$',size=17)
plt.tick_params(length=4, width=0.8, top=False, right=False, labelsize=14)
    
plt.legend(fontsize=14)
plt.show()












