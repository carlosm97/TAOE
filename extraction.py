#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 2021

@author: Martínez-Sebastián

GOAL : Extraction of the galaxy flux to get correct photometry of SN2021 aaxs.
"""

import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from astropy.time import Time
os.chdir('/home/carlos/Desktop/MSc/Segundo/TAOE/Observations')

# Firstly, we write the instrumental magnitude of our SN in different filters, different days

SN2021aaxs_21oct_R, e_21oct_R = -11.180, 0.018
SN2021aaxs_21oct_G, e_21oct_G = -10.847, 0.023
SN2021aaxs_5nov_R, e_5nov_R = -10.9114, 0.012
SN2021aaxs_5nov_G, e_5nov_G = -10.383, 0.011
SN2021aaxs_13nov_R, e_13nov_R = -14.418, 0.003
SN2021aaxs_13nov_G, e_13nov_G = -14.145, 0.002

# Zero points:
    
zp_21oct_R, ezp_21oct_R = 28.618, 0.023
zp_21oct_G, ezp_21oct_G = 28.400, 0.019
zp_5nov_R, ezp_5nov_R = 28.519, 0.021
zp_5nov_G, ezp_5nov_G = 28.287, 0.036
zp_13nov_R, ezp_13nov_R = 32.047, 0.017
zp_13nov_G, ezp_13nov_G = 32.205, 0.019

# AB magnitude: maG_ins + zp 

m21_oct_R, m21_oct_G = SN2021aaxs_21oct_R+zp_21oct_R, SN2021aaxs_21oct_G+zp_21oct_G
m5_nov_R, m5_nov_G = SN2021aaxs_5nov_R+zp_5nov_R, SN2021aaxs_5nov_G+zp_5nov_G
m13_nov_R, m13_nov_G = SN2021aaxs_13nov_R+zp_13nov_R, SN2021aaxs_13nov_G+zp_13nov_G

# errors:

dm21_oct_R, dm21_oct_G = np.sqrt(e_21oct_R**2+ezp_21oct_R**2), np.sqrt(e_21oct_G**2+ezp_21oct_G**2)
dm5_nov_R, dm5_nov_G = np.sqrt(e_5nov_R**2+ezp_5nov_R**2), np.sqrt(e_5nov_G**2+ezp_5nov_G**2)
dm13_nov_R, dm13_nov_G = np.sqrt(e_13nov_R**2+ezp_13nov_R**2), np.sqrt(e_13nov_G**2+ezp_13nov_G**2)

# Passing to mJy: 
    
flux = lambda m: 10.**(-0.4 * (m - 23.9))
    
F_R_21oct_total = flux(m21_oct_R) 
F_G_21oct_total = flux(m21_oct_G) 
F_R_5nov_total = flux(m5_nov_R) 
F_G_5nov_total = flux(m5_nov_G) 
F_R_13nov_total = flux(m13_nov_R) 
F_G_13nov_total = flux(m13_nov_G) 
del m21_oct_R,m21_oct_G,m5_nov_R,m5_nov_G,m13_nov_R,m13_nov_G

# Kron photometry of the host galaxy from PanSTARSS:
# RA, dec, radius = 128.39588031663618, 19.74125370218962, 2"

mG_Gal, dmG_Gal = 19.087, 0.091 
mR_Gal, dmR_Gal = 18.697, 0.089

# Galaxy flux: 
F_R_Gal = flux(mR_Gal) 
F_G_Gal = flux(mG_Gal) 

# SNaaxs flux:
F_R_21oct_SN = F_R_21oct_total - F_R_Gal
F_G_21oct_SN = F_G_21oct_total - F_G_Gal
F_R_5nov_SN = F_R_5nov_total - F_R_Gal
F_G_5nov_SN = F_G_5nov_total - F_G_Gal
F_R_13nov_SN = F_R_13nov_total - F_R_Gal
F_G_13nov_SN = F_G_13nov_total - F_G_Gal

del F_R_21oct_total,F_R_Gal,F_G_21oct_total,F_G_Gal,F_R_5nov_total,F_G_5nov_total,F_R_13nov_total,F_G_13nov_total

# final magnitude:
mag = lambda F: 23.9 - 2.5*np.log10(F)

m_R_21oct_SN = mag(F_R_21oct_SN)
m_G_21oct_SN = mag(F_G_21oct_SN)
m_R_5nov_SN = mag(F_R_5nov_SN)
m_G_5nov_SN = mag(F_G_5nov_SN)
m_R_13nov_SN = mag(F_R_13nov_SN)
m_G_13nov_SN = mag(F_G_13nov_SN)

# Observation MJD: http://www.csgnetwork.com/julianmodifdateconv.html
#oct21 = 59509
#nov5 = 59524
#nov13 = 59532

times = ['2021-10-22T04:59:40.856','2021-10-22T04:48:59.927',\
         '2021-11-06T04:45:13.097','2021-11-06T04:55:45.132',\
         '2021-11-14T05:02:28.8','2021-11-14T05:12:57.0']
t = Time(times, format='isot', scale='utc')
MJD_oct21_G, MJD_oct21_R, MJD_nov5_G, MJD_nov5_R, MJD_nov13_G, MJD_nov13_R =\
    t.MJD
    
# Write the data to a txt file 
try: os.remove('./SN2021aaxs_photometry.txt')
except: pass
tfile = open('SN2021aaxs_photometry.txt','a')
df = pd.DataFrame({'# date_filter': ['21_oct_R','21_oct_G','5_nov_R','5_nov_G','13_nov_R','13_nov_G'],
                   'MJD': [MJD_oct21_R,MJD_oct21_G,MJD_nov5_R,MJD_nov5_G,MJD_nov13_R,MJD_nov13_G],
                   'magnitude': [m_R_21oct_SN,m_G_21oct_SN,m_R_5nov_SN,m_G_5nov_SN,m_R_13nov_SN,m_G_13nov_SN],
                   'Error_magnitude': [dm21_oct_R,dm21_oct_G,dm5_nov_R,dm5_nov_G,dm13_nov_R,dm13_nov_G]})
tfile.write(df.to_string())
tfile.close()

# Now, we have the file to plot with TOPCAT