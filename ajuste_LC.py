#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 23 21:31:10 2022

@author: carlos
"""

import os
import pandas as pd
import numpy as np
from matplotlib import gridspec#,rcParams
import matplotlib.pylab as plt
#plt.rcParams['text.usetex'] = True
#import seaborn as sns
from astropy.time import Time
from scipy.interpolate import interp1d
os.chdir('/home/carlos/Desktop/PRUEBAS_TAOE')
# ___________________________________DATA_____________________________________
# Cargamos los datos y modelos 

#os.chdir('/home/carlos/Desktop/MSc/Segundo/TAOE/Observations/zex')
os.chdir('/home/carlos/Desktop/PRUEBAS_TAOE')
SN = 'aaxs'

# SN data 
try: data_r=pd.read_csv('SN'+SN+'_r.txt',delimiter= '\s+', index_col=False,header=0)
except: pass
data_g=pd.read_csv('SN'+SN+'_g.txt',delimiter= '\s+', index_col=False,header=0)

# Dias de observación
dia = '21oct','5nov','13nov'
#dia = ['2oct']
obs_day = {'2oct':r'${{2\,Oct}}$','21oct':r'${{21\,Oct}}$',
      '5nov':r'${{5\,Nov}}$','13nov':r'${13\,Nov}$',
      '20nov':r'${{20\,Nov}}$'}

obs_day_mjd = {'2oct':59489,'21oct':59508,'5nov':59523,'13nov':59531,'20nov':59545}
dias_mjd = [obs_day_mjd[d] for d in dia]
dias = [obs_day[d] for d in dia]



# ATLAS
data_atlas=pd.read_csv('ATLAS_'+SN+'.txt',delimiter= '\s+', index_col=False,header=0)

# ZTF
data_ztf_det=pd.read_csv('detections.csv',delimiter= ',', index_col=False,header=0)
data_ztf_nondet=pd.read_csv('non_detections.csv',delimiter= ',', index_col=False,header=0)


os.chdir("/home/carlos/Desktop/MSc/Segundo/TAOE/Observations/SN_LC")

sn2l = pd.read_csv('sn2l_lc.v1.2.dat',delimiter= '\s+',index_col=False,header=0)
sn2p = pd.read_csv('sn2p_lc.v1.2.dat',delimiter= '\s+',index_col=False,header=0)
sn2n = pd.read_csv('sn2n_lc.v2.1.dat',delimiter= '\s+',index_col=False,header=0)
    
'''
sn2l['m'] /= np.nanmax(sn2l['m'])    
sn2p['m'] /= np.nanmax(sn2p['m'])    
sn2n['m'] /= np.nanmax(sn2n['m'])    
'''
#%%
# Buscamos un primer ajuste "automático" a la LC.
plt.close('all')
#fig=plt.figure(figsize=(10,7))  
spec = gridspec.GridSpec(ncols=1, nrows=2,height_ratios=[4,1])
 
ax0 = fig.add_subplot(spec[0,0])
plt.gca().invert_yaxis()
plt.grid()
# ZTF rojo
plt.errorbar(data_ztf_det['mjd'].mask(data_ztf_det['fid']==1),data_ztf_det['magpsf'].mask(data_ztf_det['fid']==1),\
             yerr=data_ztf_det['sigmapsf'].mask(data_ztf_det['fid']==1),fmt='o',ls='none',markerfacecolor='white',color='red',label='ZTF r')
# ZTF verde
plt.errorbar(data_ztf_det['mjd'].mask(data_ztf_det['fid']==2),data_ztf_det['magpsf'].mask(data_ztf_det['fid']==2),\
             yerr=data_ztf_det['sigmapsf'].mask(data_ztf_det['fid']==2),fmt='o',ls='none',markerfacecolor='white',color='lime',label='ZTF g')

# TAOE
try: plt.errorbar(data_r['mjd'],data_r['m'],yerr=data_r['um'],fmt='o',ls='none',color='red',label='TAOE r')
except: pass
plt.errorbar(data_g['mjd'],data_g['m'],yerr=data_g['um'],fmt='o',ls='none',color='green',label='TAOE g')
    
ref_magnitude = min(np.nanmin(data_ztf_det['magpsf'].mask(data_ztf_det['fid']==1)),np.nanmin(data_ztf_det['magpsf'].mask(data_ztf_det['fid']==2)))
ref_MJD = data_ztf_det['mjd'][data_ztf_det['magpsf']==ref_magnitude]

plt.plot(sn2l['d']+float(ref_MJD),sn2l['m']+ref_magnitude,label='SN IIl V')
plt.plot(sn2p['d']+float(ref_MJD),sn2p['m']+ref_magnitude,label='SN IIp V')
plt.plot(sn2n['d']+float(ref_MJD)-20,sn2n['m']+ref_magnitude,label='SN IIn V')
plt.ylabel('${m}$',size=17),plt.xlabel('${MJD}$',size=17)
plt.tick_params(length=4, width=0.8, top=False, right=False, labelsize=14)
plt.legend(fontsize=14)

plt.ylim(np.nanmax(data_ztf_det['magpsf'])+0.5,np.nanmin(data_ztf_det['magpsf'])-0.5)
plt.xlim(np.nanmin(data_ztf_det['mjd'])-10,np.nanmax(data_ztf_det['mjd'])+10)
















    
    
    
    
    
    
