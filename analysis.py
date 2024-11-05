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
'''
read_files = glob.glob("*.txt")

with open("result.txt", "wb") as outfile:
    for f in read_files:
        with open(f, "rb") as infile:
            outfile.write(infile.read())

mpl.rc('font', family='serif', serif='Linguistics Pro')
mpl.rc('text', usetex=False)
mpl.rc('mathtext', fontset='custom',
       rm='Linguistics Pro',
       it='Linguistics Pro:italic',
       bf='Linguistics Pro:bold')
'''
# Define some test data which is close to Gaussian
column_names = pd.read_csv('column_names.txt',sep='	', float_precision='high')
data_whole = pd.read_csv('20140807165108.txt',sep='	', float_precision='high')
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

'''
n=10 # number of points to be checked before and after 
# Find local peaks
df['max'] = df.iloc[argrelextrema(df.rol.values, np.greater_equal, order=n)[0]]['rol']
peaks, _ = find_peaks(df.rol.values, height=0)
#plt.scatter(df.latitude, df['max'], c='g')
#ax.set_ylim(15000,25000)
#data_whole.rename(columns=lambda s: s*3)
#data_whole.rename(columns={0:542})
#data=data_whole.loc[:, 1:5]

data_whole[['tof',"sweep", "counts"]] = pd.DataFrame([ x.split() for x in data_whole["[DATA]"].tolist() ])
data_whole = data_whole.drop('[DATA]', axis=1).astype(float)

z_cut = 5              # z_cut = number_of_max_ions
tof_window = [8000,14000]        # tof_window = [tof_min, tof_max]
#sweep_windows = [[]]   # sweep_windows = [[sweep_min_1, sweep_max_1], [sweep_min_2, sweep_max_2], ...]

# z-class cut
data_whole = data_whole[data_whole.sweep.isin(data_whole.sweep.value_counts()[data_whole.sweep.value_counts() >= z_cut].index)]
# tof-window cut
if tof_window != []:
    data_whole = data_whole[(data_whole.tof >= tof_window[0]) & (data_whole.tof <= tof_window[1])]
# sweep-windows cut
  
fig, ((ax_x, blank),(ax_0, ax_y)) = plt.subplots(2,2,sharex='col',sharey='row', figsize=(9,9),
                                                 gridspec_kw={'height_ratios':[1,4],
                                                             'width_ratios':[4,1],
                                                             'hspace': 0.05,
                                                             'wspace':0.05})
blank.remove()

data_whole.plot(x='tof', y='sweep', style='o', alpha=0.15, ms=2,ax=ax_0, label='unbinned data')
ax_0.set_xlabel('time of flight (ns)', fontsize=24, fontweight='bold')
ax_0.set_ylabel('sweep', fontsize=24, fontweight='bold')

x_proj = data_whole.tof.value_counts(bins=500).sort_index()
y_proj = data_whole.sweep.value_counts(bins=500).sort_index()
bin_centers=x_proj.index.mid
np.savetxt("hist.txt", zip(x_proj.index.left, x_proj.index.right, x_proj))
bin_data = pd.read_csv('hist.txt', sep=" ", float_precision='high')
bin_data.columns = ['left',"right", "c1"]

ax_x.bar(x_proj.index.mid, x_proj, width = 15)
#ax_x.set_yscale('log')
ax_y.plot(y_proj, y_proj.index.mid)

ax_0.set_xlabel('time of flight (ns)', fontsize=24, fontweight='bold')
ax_0.set_ylabel('sweep', fontsize=24, fontweight='bold')
# ax_0.ticklabel_format(useOffset=False, style='plain')
ax_x.set_ylabel('# / 4 ns', fontsize=24, fontweight='bold')
ax_y.xaxis.set_ticks_position('top')
ax_y.xaxis.set_label_position('top')

peak_threshold = 0.9
peak_min_distance = 70
s=3

x_proj_peak_ind = peakutils.indexes(bin_data['c1'], thres=peak_threshold, min_dist=peak_min_distance)
peak = pd.DataFrame(bin_data, columns=bin_data.columns, index=x_proj_peak_ind)
x_proj_peak_inter = peakutils.interpolate(bin_data.index, bin_data['c1'], ind=x_proj_peak_ind)
peak['midele']=(peak.left+peak.right)/2
a_list = []
b_list = []
c_list = []
d_list = []
for peak_pos in peak['midele']:
    ax_x.axvline(peak_pos, c='g')
    ax_x.axvline(peak_pos-300, c='g', alpha=0.3)
    ax_x.axvline(peak_pos+300, c='g', alpha=0.3)
    ax_0.axvline(peak_pos, c='g', label='peak position' if peak_pos == x_proj_peak_ind[0] else "")
    ax_0.axvline(peak_pos-300, c='g', alpha=0.3)
    ax_0.axvline(peak_pos+300, c='g', alpha=0.3)

for (peak_left, peak_right) in zip(peak['left'], peak['right']):
    for (x,y,z) in zip(data_whole['tof'],data_whole['sweep'],data_whole['counts']):
        if peak_left-300<=x<=peak_right+300:
            a=x
            b=y
            c=z
            d=(peak_left+peak_right)/2
            a_list.append(a)
            b_list.append(b)
            c_list.append(c)
            d_list.append(d)
important = pd.DataFrame({'tof': a_list, 'sweeps': b_list, 'counts': c_list, 'peak': d_list})
del a_list, b_list, c_list, d_list
df_list_1 = [g for _, g in important.groupby(['peak'])]
nr=10
std_value_0=df_list_1[0]
std_value_0=std_value_0['tof'].rolling(nr).mean().std()

std_value_1=df_list_1[1]
std_value_1=std_value_1['tof'].rolling(nr).mean().std()
std_value_2=df_list_1[2]
std_value_2=std_value_2['tof'].rolling(nr).mean().std()

data=[std_value_0]#,std_value_1,std_value_2]
std_value = pd.DataFrame(data, columns = ['std'])
a_list = []
b_list = []
c_list = []
d_list = []
for (std_value, peak_middle) in zip(std_value['std'], peak['midele']):
    for (x,y,z) in zip(data_whole['tof'],data_whole['sweep'],data_whole['counts']):
        if peak_middle-s*std_value<=x<=peak_middle+s*std_value:
            a=x
            b=y
            c=z
            d=peak_middle
            a_list.append(a)
            b_list.append(b)
            c_list.append(c)
            d_list.append(d)
filtered = pd.DataFrame({'tof': a_list, 'sweeps': b_list, 'counts': c_list, 'peak': d_list})
del a_list, b_list, c_list, d_list
df_list_2 = [g for _, g in filtered.groupby(['peak'])]
x_proj_filtered = filtered.tof.value_counts(bins=500).sort_index()
ax_x.bar(x_proj_filtered.index.mid, x_proj_filtered, width = 15,color='red')
#ax_x.set_yscale('log')
filtered.plot(x='tof', y='sweeps', style='o', alpha=0.15, ms=2,ax=ax_0, label='unbinned data',c='red')
tof_value_0=df_list_1[0]
tof_value_0=tof_value_0['tof']
box_plot=[tof_value_0,np.repeat(peak.midele,len(tof_value_0)),tof_value_0+s*std_value,tof_value_0-s*std_value]
#box_plot=[tof_value_0,np.repeat(peak.midele,len(tof_value_0)),np.repeat(peak.midele+3*std_value,len(tof_value_0)),np.repeat(peak.midele-1*std_value,len(tof_value_0))]
box_plot[1].index=box_plot[2].index=box_plot[3].index=box_plot[0].index
box_data = np.concatenate((box_plot[0], box_plot[1], box_plot[2], box_plot[3]))
#box_plot=np.concatenate()
fig5, ax5 = plt.subplots()
ax5.set_title('Horizontal Boxes')
#red_square = dict(markerfacecolor='r', marker='s')
ax5.boxplot(box_data, vert=False)


t = mle.var('t', observed=True, vector=True)
mean_t = mle.var('mean_t')
sigma_t = mle.var('sigma_t')
# defining model
model = mle.Normal(t, mean_t, sigma_t)
# defining fit ranges
fit_ranges = [40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 110, 120, 130, 140, 160, 180]
results_df_dict = {}
for i in range(len(x_proj_peak_ind)):
    results_df_dict['peak_{}'.format(i+1)] = []
    for fit_range in fit_ranges:
        # applying fit range by cutting the data frame
        fit_df = cut_df.tof[(cut_df.tof < x_proj_peak_ind[i] + fit_range) &
                            (cut_df.tof > x_proj_peak_ind[i] - fit_range)]
        # fit for all peaks and different ranges
        r = model.fit({'t': fit_df.to_numpy()},
                     {'mean_t': fit_df.mean(),
                      'sigma_t': 60.0})
        results_df_dict['peak_{}'.format(i+1)].append([r['x']['mean_t'], np.sqrt(r['hess_inv'][0][0]) , r['x']['sigma_t'], np.sqrt(r['hess_inv'][1][1])])

# save fit results in dataframe
for key in results_df_dict.keys():
    results_df_dict[key] = pd.DataFrame(results_df_dict[key], columns=['mean_t', 'mean_t_unc', 'sigma', 'sigma_unc'])
results_df = pd.concat(results_df_dict, axis=1).set_index(pd.Series(fit_ranges))

results_df.to_csv('data/{}-fits-summary.csv'.format(file))


for (sweeps,tof,counts,peak) in zip(importan['sweeps'],importan['tof'],importan['counts'],importan['peak']):
    if (sweeps==sweeps+1) and (peak==peak+1):
        a=
        a_list.append(a)
        b_list.append(b)
        c_list.append(c)
        d_list.append(d)
        
    sweep_outs = []
    sweep_ins = []
        # find positions where rolling avg crosses +/-2 sigma of the main peak outwards to exclude these files
    for sweep_out in (data_whole.sweep[(data_whole.tof < peak_pos+300) & (data_whole.tof > peak_pos-300)]
                              [((rolling_ser.fillna(value=peak_pos) - (peak_pos+2*rolling_std) > 0) &
                                 (rolling_ser.fillna(value=peak_pos).shift() - (peak_pos+2*rolling_std) < 0)) |
                               ((rolling_ser.fillna(value=peak_pos) - (peak_pos-2*rolling_std) < 0) &
                                (rolling_ser.fillna(value=peak_pos).shift() - (peak_pos-2*rolling_std) > 0))]
                          ):
        sweep_outs.append(sweep_out)
        # find positions where rolling avg crosses +/-2 sigma of the main peak inwards to exclude these files
    for sweep_in in (data_whole.sweep[(data_whole.tof < peak_pos+300) & (data_whole.tof > peak_pos-300)]
                              [((rolling_ser.fillna(value=peak_pos) - (peak_pos+2*rolling_std) < 0) &
                                (rolling_ser.fillna(value=peak_pos).shift() - (peak_pos+2*rolling_std) > 0)) |
                               ((rolling_ser.fillna(value=peak_pos) - (peak_pos-2*rolling_std) > 0) &
                                (rolling_ser.fillna(value=peak_pos).shift() - (peak_pos-2*rolling_std) < 0))]
                          ):
        sweep_ins.append(sweep_in)
        # plot areas where rolling average is outside the +/- 2 sigma band
    for i in range(len(sweep_outs)):
                ax_0.axhspan(sweep_outs[i], sweep_ins[i], color='k', alpha=0.2,
                         label='statistical cut' if i == 0 else "")
                ax_y.axhspan(sweep_outs[i], sweep_ins[i], color='k', alpha=0.2)


   
    if peak_pos in x_proj[x_proj_peak_ind].idxmax():
        rolling_nr = 350
        # calculates the standard deviation of the rolling average
        rolling_std = data_whole.tof[(data_whole.tof < peak_pos+300) &
                    (data_whole.tof > peak_pos-300)].rolling(rolling_nr).mean().std()
        # calculates the rolling average and shifts it by half of rolling_nr for proper display
        rolling_ser = data_whole.tof[(data_whole.tof < peak_pos+300) &
                    (data_whole.tof > peak_pos-300)].rolling(rolling_nr).mean().shift(-int(rolling_nr/2))
        # plotting the rolling average       
        ax_0.plot(rolling_ser, data_whole.sweep[(data_whole.tof < peak_pos+300) &
                    (data_whole.tof > peak_pos-300)], '-', lw=0.5, c='r',
                  label='rolling avg. (n={})'.format(rolling_nr))
        ax_0.axvline(peak_pos+2*rolling_std, c='r', alpha=0.5)
        ax_0.axvline(peak_pos-2*rolling_std, c='r', alpha=0.5)
        sweep_outs = []
        sweep_ins = []
        # find positions where rolling avg crosses +/-2 sigma of the main peak outwards to exclude these files
        for sweep_out in (data_whole.sweep[(data_whole.tof < peak_pos+300) & (data_whole.tof > peak_pos-300)]
                              [((rolling_ser.fillna(value=peak_pos) - (peak_pos+2*rolling_std) > 0) &
                                 (rolling_ser.fillna(value=peak_pos).shift() - (peak_pos+2*rolling_std) < 0)) |
                               ((rolling_ser.fillna(value=peak_pos) - (peak_pos-2*rolling_std) < 0) &
                                (rolling_ser.fillna(value=peak_pos).shift() - (peak_pos-2*rolling_std) > 0))]
                          ):
            sweep_outs.append(sweep_out)
        # find positions where rolling avg crosses +/-2 sigma of the main peak inwards to exclude these files
        for sweep_in in (data_whole.sweep[(data_whole.tof < peak_pos+300) & (data_whole.tof > peak_pos-300)]
                              [((rolling_ser.fillna(value=peak_pos) - (peak_pos+2*rolling_std) < 0) &
                                (rolling_ser.fillna(value=peak_pos).shift() - (peak_pos+2*rolling_std) > 0)) |
                               ((rolling_ser.fillna(value=peak_pos) - (peak_pos-2*rolling_std) > 0) &
                                (rolling_ser.fillna(value=peak_pos).shift() - (peak_pos-2*rolling_std) < 0))]
                          ):
            sweep_ins.append(sweep_in)
        # plot areas where rolling average is outside the +/- 2 sigma band
        for i in range(len(sweep_outs)):
            ax_0.axhspan(sweep_outs[i], sweep_ins[i], color='k', alpha=0.2,
                         label='statistical cut' if i == 0 else "")
            ax_y.axhspan(sweep_outs[i], sweep_ins[i], color='k', alpha=0.2)


data = pd.read_csv('Try_1.txt', sep=" ", header=None, float_precision='high')
data.columns = ['space1',"x1", "counts1"]
bin_centers=a1=data['x1'].astype(float)#*1e-6
b1=data['counts1'].astype(float)

#plt.figure(figsize=[10,8])
plt.bar(a1, b1, width = 0.9, color='#0504aa',alpha=0.7)
A=100
mu= 16856447
sigma=21.27
B=numpy.exp(-(a1-mu)**2/(2.*sigma**2))
gaus_check=A*numpy.exp(-(bin_centers-mu)**2/(2.*sigma**2))
print(gaus_check)



fig, ax = plt.subplots(1)
#plt.xlim(min(bin_centers), max(bin_centers))

ax.ticklabel_format( axis='x',style='plain',useOffset=False)
ax.set_ylabel('Counts/0.8 ns',fontsize=16)
ax.set_xlabel('Time of flight (ms)',fontsize=16)
ax.tick_params(left=True, right=True, top=True, bottom=True, labelleft='on', labelright='on', labelsize=12)
ax.tick_params(which='minor', right=True)
ax.set_yscale('log')
#ax.set_ylim(0,35)
ax.set_xlim(min(a),max(a))
ax.xaxis.set_ticks(np.arange(min(a),max(a), (max(a)-min(a))/4.5))


def gauss(x, *p):
    A, mu, sigma = p
    return A*numpy.exp(-(x-mu)**2/(2.*sigma**2))
p0 = [1000000*1e-6, 20575385*1e-6, 1.6*1e-6]
coeff, var_matrix = curve_fit(gauss, bin_centers, b1, p0=p0)
hist_gaus_fit = gauss(bin_centers, *coeff)
ax.plot(bin_centers, hist_gaus_fit, color='red',)
print(coeff[1])


def d_gauss(x, *p):
    A1, mu1, sigma1, A2, mu2, sigma2 = p
    return A1*numpy.exp(-(x-mu1)**2/(2.*sigma1**2))+A2*numpy.exp(-(x-mu2)**2/(2.*sigma2**2))
p0 = [3000000*1e-6, 20575280.2*1e-6, 3.6*1e-6,1000000*1e-5, 20575387.6*1e-6, 1.7*1e-6]
coeff, var_matrix = curve_fit(d_gauss, bin_centers, b1, p0=p0)
hist_d_gaus_fit = d_gauss(bin_centers, *coeff)
ax.plot(bin_centers, hist_d_gaus_fit, color='red',)
print(coeff)

p0 = [1000000*1e-6, 16856079*1e-6, 19.98*1e-6]
coeff_ti, var_matrix_ti = curve_fit(gauss, bin_centers, b1, p0=p0)
hist_gaus_fit_ti = gauss(bin_centers, *coeff_ti)
ax.plot(bin_centers, hist_gaus_fit_ti, color='blue',)
#plt.axvline(x=coeff[1],color='black',linestyle='--')


def egh(x, *p1):
    A1, mu1, sigma1, tau1 = p1
    return A1*(tau1/2)*numpy.exp(tau1*(2*mu1+tau1*sigma1**2-2*x)/2)*erfc((mu1+tau1*sigma1**2-x)/(sigma1*np.sqrt(2)))
lamda=1/(1.77*1e-6)
p0 = [coeff[0], 16856447*1e-6, 21.79*1e-6, lamda]
coeff1, var_matrix1 = curve_fit(egh, bin_centers, b1, p0=p0)
hist_egh_fit = egh(bin_centers, *coeff1)
ax.plot(bin_centers, hist_egh_fit, color='blue',)
print(coeff1[1])
lamda_Ti=1/(1.87*1e-6)
p0 = [coeff[0], 16856078*1e-6, 19.92*1e-6, lamda_Ti]
coeff_Ti, var_matrix_Ti = curve_fit(egh, bin_centers, b1, p0=p0)
hist_egh_fit_Ti = egh(bin_centers, *coeff_Ti)
ax.plot(bin_centers, hist_egh_fit_Ti, color='red',)
#plt.axvline(x=coeff1[1],color='red',linestyle='--')

#ax.legend(('Gaussian'),loc='upper left')
ax.bar(a, b, width = 0.000005, color='k',alpha=0.7)
fig.text(0.26,0.779, '$^{73}$Ge$^+$ \n $^{73}$As$^+$ ', ha='center', va='center', fontsize=16)
fig.text(0.335,0.655, '$^{73}$Ga$^+$ ', ha='center', va='center', fontsize=16)
fig.text(0.435,0.82, '$^{73}$Se$^+$ ', ha='center', va='center', fontsize=16)
fig.text(0.72,0.52, '$^{61}$Ni$^{12}$C$^+$ ', ha='center', va='center', fontsize=16)
fig.text(0.82,0.67, '$^{73}$Br$^+$ ', ha='center', va='center', fontsize=16)
ax.legend(('fit','data'),loc='upper right')

plt.savefig("fit.pdf",  bbox_inches='tight')


'''
#plt.savefig("1.pdf",  bbox_inches='tight')
plt.show()
