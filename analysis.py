# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 12:29:57 2019

@author: ikulikov
"""

from scipy import stats
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import numpy as np
import math
from scipy.signal import argrelextrema
from scipy.signal import find_peaks
#import pickle, os
#import math
#import mle
import glob
import peakutils

column_names = pd.read_csv('column_names.txt',sep='	', float_precision='high')
data_whole = pd.read_csv('20141013040046.txt',sep='	', float_precision='high')
data_whole.columns=column_names.columns
data=pd.concat([data_whole.iloc[:,0:6], data_whole.iloc[:,7:262].sum(axis=1)], axis=1)
data=pd.concat([data, data_whole.iloc[:,519:531]], axis=1)
data.columns=['year','month','day','hour','minute','second','intensity','latitude','longitude','521','522','523','524','525','526','527','528','529','L']
a=4.6*10**(24)/data_whole['6']**8.3
data=pd.concat([data, a*data['intensity']/256], axis=1)
fig, ax = plt.subplots()
ax.set_title('Data')
nr_1=2
nr_2=4
data['std']=data[0].rolling(nr_1).mean().std()
data['rol']=data[0].rolling(nr_1).mean()
ax.plot(data['latitude'],data['rol'],color='red',alpha=0.5 )
ax.set_yscale('log')
s=6#number of standart deviation from mean
df = data.copy(deep = True)
df['6']=data_whole['6']
df = df[data[0] < s*df['std']]
df['rol']=df[0].rolling(nr_2).mean()
ax.plot(df['latitude'],df['rol'],color='blue')

threshold=0.1
minimum_distanse=28
cb = np.array(data.rol.dropna())
peaks = peakutils.indexes(cb, thres=threshold, min_dist=minimum_distanse)
peaks= pd.DataFrame(data, index=peaks)
ax.scatter(peaks.latitude,peaks.rol,c='orange',marker='o',alpha=1)

plt.show()
