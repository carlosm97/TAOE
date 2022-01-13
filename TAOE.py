# -*- coding: utf-8 -*-
"""
Created on Sun Apr 25 14:13:22 2021

@author: Mar Pérez Sar

GOAL : Graficar ATLAS + ZTF y curva de luz. 
"""


import os
import pandas as pd
import numpy as np
from matplotlib import rcParams,gridspec
import matplotlib.pylab as plt
#plt.rcParams['text.usetex'] = True
#import seaborn as sns
from astropy.time import Time
from scipy.interpolate import interp1d
os.chdir('/home/carlos/Desktop/PRUEBAS_TAOE')
# ___________________________________DATA_____________________________________

# Dias de observación
dia = '21oct','5nov','13nov'

obs_day = {'2oct':r'${{2\,Oct}}$','21oct':r'${{21\,Oct}}$',
      '5nov':r'${{5\,Nov}}$','13nov':r'${13\,Nov}$',
      '20nov':r'${{20\,Nov}}$'}

obs_day_mjd = {'2oct':59489,'21oct':59508,'5nov':59523,'13nov':59531,'20nov':59545}
dias_mjd = [obs_day_mjd[d] for d in dia]
dias = [obs_day[d] for d in dia]

SN = 'aaxs'

# ATLAS
data_atlas=pd.read_csv('ATLAS_'+SN+'.txt',delimiter= '\s+', index_col=False,header=0)

# ZTF
data_ztf_det=pd.read_csv('detections.csv',delimiter= ',', index_col=False,header=0)
data_ztf_nondet=pd.read_csv('non_detections.csv',delimiter= ',', index_col=False,header=0)

# SN data 

data_r=pd.read_csv('SN'+SN+'_r.txt',delimiter= '\s+', index_col=False,header=0)
data_g=pd.read_csv('SN'+SN+'_g.txt',delimiter= '\s+', index_col=False,header=0)

#__________________________________________________________________________________

# PRIMEIRA FIGURA. COMPARACIÓN ATLAS VS ZTF
fig=plt.figure(figsize=(14,7))  
# plt.suptitle('$\displaystyle{\mathrm{{SNxfu}}}$',size=20)
spec = gridspec.GridSpec(ncols=2, nrows=2,width_ratios=[1,2]) 

# ATLAS  
ax0 = fig.add_subplot(spec[0,0])
plt.title(r'${{ATLAS}}$',size=17,y=.05, x=0.85)
# Filtro naranja
plt.errorbar(data_atlas['###MJD'].mask(data_atlas['F']=='o'),data_atlas['m'].mask(data_atlas['F']=='o'),yerr=data_atlas['dm'],fmt='o',markersize=5,ls='none',color='cyan',label='Cyan')
# Filtro cyan
plt.errorbar(data_atlas['###MJD'].mask(data_atlas['F']=='c'),data_atlas['m'].mask(data_atlas['F']=='c'),yerr=data_atlas['dm'],fmt='o',markersize=5,ls='none',color='darkorange',label='Orange')
plt.grid(which='major', axis='both',alpha=0.3, linestyle='-')
plt.setp(ax0.get_xticklabels(), visible=False)
plt.ylabel('${m}$',size=17)
plt.tick_params(length=4, width=0.8, top=False, right=False, labelsize=14)
#plt.ylim(int(2*min(np.nanmedian(data_atlas['m'].mask(data_atlas['F']=='o'))-2*np.std(data_atlas['m'].mask(data_atlas['F']=='o')),np.nanmedian(data_atlas['m'].mask(data_atlas['F']=='c'))-2*np.std(data_atlas['m'].mask(data_atlas['F']=='c'))))/2,
#         int(2*max(np.nanmedian(data_atlas['m'].mask(data_atlas['F']=='o'))+2*np.std(data_atlas['m'].mask(data_atlas['F']=='o')),np.nanmedian(data_atlas['m'].mask(data_atlas['F']=='c'))+2*np.std(data_atlas['m'].mask(data_atlas['F']=='c')))+1)/2)
plt.ylim(19.5,17)
plt.gca().invert_yaxis()

# ZTF
ax1 = fig.add_subplot(spec[1,0],sharex=ax0)
plt.title('$ZTF$',size=17,y=.05, x=0.86)
# Filtro verde
plt.errorbar(data_ztf_det['mjd'].mask(data_ztf_det['fid']==2),data_ztf_det['magpsf'].mask(data_ztf_det['fid']==2),yerr=data_ztf_det['sigmapsf'],fmt='o',markersize=5,ls='none',color='lime',label='Green')
plt.plot(data_ztf_nondet['mjd'].mask(data_ztf_nondet['fid']==2),data_ztf_nondet['diffmaglim'].mask(data_ztf_nondet['fid']==2),'v',color='lime',alpha=0.5)
# Filtro rojo
plt.errorbar(data_ztf_det['mjd'].mask(data_ztf_det['fid']==1),data_ztf_det['magpsf'].mask(data_ztf_det['fid']==1),yerr=data_ztf_det['sigmapsf'],fmt='o',markersize=5,ls='none',color='red',label='Red')
plt.plot(data_ztf_nondet['mjd'].mask(data_ztf_nondet['fid']==1),data_ztf_nondet['diffmaglim'].mask(data_ztf_nondet['fid']==1),'v',color='red',alpha=0.5)
plt.grid(which='major', axis='both',alpha=0.3, linestyle='-')
plt.xlabel(r'${MJD}$',size=17)
plt.ylabel(r'${m}$',size=17)
#plt.ylim(18,22)
plt.tick_params(length=4, width=0.8, top=False, right=False, labelsize=14)
plt.gca().invert_yaxis()
plt.subplots_adjust(hspace=0.09)


# COMBINACION
ax2 = fig.add_subplot(spec[:,1])

# ATLAS
# Filtro naranja
plt.plot(data_atlas['###MJD'].mask(data_atlas['F']=='o'),data_atlas['m'].mask(data_atlas['F']=='o'),'o',markersize=5,ls='none',color='cyan')
# Filtro cyan
plt.plot(data_atlas['###MJD'].mask(data_atlas['F']=='c'),data_atlas['m'].mask(data_atlas['F']=='c'),'o',markersize=5,ls='none',color='darkorange')


# ZTF
# Filtro verde
plt.errorbar(data_ztf_det['mjd'].mask(data_ztf_det['fid']==2),data_ztf_det['magpsf'].mask(data_ztf_det['fid']==2),yerr=data_ztf_det['sigmapsf'],fmt='o',ls='none',color='lime')
#plt.plot(data_ztf_nondet['mjd'].mask(data_ztf_nondet['fid']==2),data_ztf_nondet['diffmaglim'].mask(data_ztf_nondet['fid']==2),'v',color='lime',alpha=0.5)
# Filtro rojo
plt.errorbar(data_ztf_det['mjd'].mask(data_ztf_det['fid']==1),data_ztf_det['magpsf'].mask(data_ztf_det['fid']==1),yerr=data_ztf_det['sigmapsf'],fmt='o',ls='none',color='red')
#plt.plot(data_ztf_nondet['mjd'].mask(data_ztf_nondet['fid']==1),data_ztf_nondet['diffmaglim'].mask(data_ztf_nondet['fid']==1),'v',color='red',alpha=0.5)
#plt.xlim(59421,59538+5)
plt.ylim(int(2*min(np.nanmedian(data_ztf_det['magpsf'].mask(data_ztf_det['fid']==2))-3*np.std(data_ztf_det['magpsf'].mask(data_ztf_det['fid']==2)),np.nanmedian(data_ztf_det['magpsf'].mask(data_ztf_det['fid']==1))-3*np.std(data_ztf_det['magpsf'].mask(data_ztf_det['fid']==1))))/2,
         int(2*max(np.nanmedian(data_ztf_det['magpsf'].mask(data_ztf_det['fid']==2))+3*np.std(data_ztf_det['magpsf'].mask(data_ztf_det['fid']==2)),np.nanmedian(data_ztf_det['magpsf'].mask(data_ztf_det['fid']==1))+3*np.std(data_ztf_det['magpsf'].mask(data_ztf_det['fid']==1)))+1)/2)
plt.gca().invert_yaxis()



# Combinacion
plt.grid(which='major', axis='both',alpha=0.3, linestyle='-')
plt.xlabel(r'${MJD}$',size=17)
plt.ylabel(r'${m}$',size=17)
plt.tick_params(length=4, width=0.8, top=False, right=False, labelsize=14)

for d in range(len(dias_mjd)):
    plt.axvline(x=dias_mjd[d], ymin=0, ymax=1,linestyle='-.',color='black',alpha=0.6)
    plt.text(dias_mjd[d]-3,16.95,dias[d], color='black',size=14,rotation=45)
plt.subplots_adjust(wspace=0.3)


#%%
# SEGUNDA GRÁFICA. CURVA DE LUZ CON LA EVOLUCIÓN DEL COLOR ABAJO__________________

fig=plt.figure(figsize=(10,7))  
spec = gridspec.GridSpec(ncols=1, nrows=2,height_ratios=[4,1])
 
ax0 = fig.add_subplot(spec[0,0])

# ZTF rojo
plt.errorbar(data_ztf_det['mjd'].mask(data_ztf_det['fid']==1),data_ztf_det['magpsf'].mask(data_ztf_det['fid']==1),yerr=data_ztf_det['sigmapsf'].mask(data_ztf_det['fid']==1),fmt='o',ls='none',markerfacecolor='white',color='red')
# ZTF verde
plt.errorbar(data_ztf_det['mjd'].mask(data_ztf_det['fid']==2),data_ztf_det['magpsf'].mask(data_ztf_det['fid']==2),yerr=data_ztf_det['sigmapsf'].mask(data_ztf_det['fid']==2),fmt='o',ls='none',markerfacecolor='white',color='lime')


# XFU nuestra
plt.errorbar(data_r['mjd'],data_r['m'],yerr=data_r['um'],fmt='o',ls='none',color='red')
plt.errorbar(data_g['mjd'],data_g['m'],yerr=data_g['um'],fmt='o',ls='none',color='green')


plt.ylim(int(2*min(np.nanmedian(data_r['m'])-3*np.std(data_r['m']),np.nanmedian(data_g['m'])-3*np.std(data_g['m'])))/2,
         int(2*max(np.nanmedian(data_r['m'])+3*np.std(data_r['m']),np.nanmedian(data_g['m'])+3*np.std(data_g['m']))+1)/2)
plt.gca().invert_yaxis()
plt.grid(which='major', axis='both',alpha=0.3, linestyle='-')
plt.setp(ax0.get_xticklabels(), visible=False)
plt.ylabel('${m}$',size=17)
plt.tick_params(length=4, width=0.8, top=False, right=False, labelsize=14)
for d in range(len(dias_mjd)):
    plt.axvline(x=dias_mjd[d], ymin=0, ymax=1,linestyle='-.',color='black',alpha=0.6)
    plt.text(dias_mjd[d]-3,16.95, dias[d], color='black',size=14,rotation=45)


ax0 = fig.add_subplot(spec[1,0],sharex=ax0)

r=interp1d(data_ztf_det['mjd'].mask(data_ztf_det['fid']==1),data_ztf_det['magpsf'].mask(data_ztf_det['fid']==1))
g=interp1d(data_ztf_det['mjd'].mask(data_ztf_det['fid']==2),data_ztf_det['magpsf'].mask(data_ztf_det['fid']==2))
err=np.sqrt(data_ztf_det['sigmapsf'].mask(data_ztf_det['fid']==1)**2+data_ztf_det['sigmapsf'].mask(data_ztf_det['fid']==2)**2)
x=data_ztf_det['mjd'].mask(data_ztf_det['fid']==2)

plt.errorbar(x,g(x)-r(x),yerr=data_ztf_det['sigmapsf'],fmt='o',ls='none',markerfacecolor='white',color='black')

plt.grid(which='major', axis='both',alpha=0.3, linestyle='-')
plt.xlabel('${MJD}}$',size=17)
plt.ylabel('${(g-r)_{ZTF}}$',size=17)
plt.tick_params(length=4, width=0.8, top=False, right=False, labelsize=14)
for d in dias_mjd:
    plt.axvline(x=d, ymin=0, ymax=1,linestyle='-.',color='black',alpha=0.6)
plt.subplots_adjust(hspace=0.09)