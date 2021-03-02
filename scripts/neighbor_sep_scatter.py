from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import pandas as pd
import matplotlib.lines as mlines
import os.path
import scipy.stats


loc = '/Users/josh/projects/intro/sources/'
res = 150
sources = pd.read_csv('sources.csv')
sources = sources['circinus'].values.tolist()
all_first, all_second, all_third = [], [], []
for i in range(len(sources)):
    fp = loc+str(sources[i])+'/'+str(sources[i])+'_'+str(res)+'pc_cloud_stats.csv'
    if os.path.isfile(fp):
        cat = pd.read_csv(fp)

        #ALL GALAXIES
        all_first = all_first + cat['min_dist'].tolist()
        all_second = all_second + cat['min_dist2nd'].tolist()
        all_third = all_third + cat['min_dist3rd'].tolist()



num_bins = 15
bins = scipy.stats.binned_statistic(all_first, all_second, statistic='median', bins=num_bins, range=(0,1000))
medians = bins[0]
xvals = (bins[1][1:] + bins[1][:-1]) / 2

bins2 = scipy.stats.binned_statistic(all_first, all_second, statistic='median', bins=num_bins, range=(0,1000))
medians2 = bins2[0]
xvals2 = (bins2[1][1:] + bins2[1][:-1]) / 2
#Get percentiles
err1 = np.zeros((2,num_bins))
err2 = np.zeros((2,num_bins))
for i in range(num_bins):
        err1[:,i] = [np.percentile(np.take(all_second, np.where(bins[2]==(i+1))),16),np.percentile(np.take(all_second, np.where(bins[2]==(i+1))),84)]
        err2[:,i] = [np.percentile(np.take(all_third, np.where(bins[2]==(i+1))),16),np.percentile(np.take(all_third, np.where(bins[2]==(i+1))),84)]
err1 = np.abs(err1 - medians)
err2 = np.abs(err2 - medians2)
fig, axes = plt.subplots(1, 2, figsize=(10,5))
axes[0].scatter(all_first, all_second, c='gray', s=3, alpha=0.2, label=r'$\rho$ = '+ str(np.around(scipy.stats.pearsonr(all_first, all_second)[0], 2)))
axes[0].set_xlim((0,1000))
axes[0].set_ylim((0,2000))
axes[0].scatter(xvals, medians, c='b')
axes[0].errorbar(xvals, medians, c='b', capsize=4, ls='none', xerr=None, yerr=err1)
axes[0].set_title(str(res)+'pc Resolution')
axes[0].set_xlabel('Distance to 1st Nearest Neighbor (pc)')
axes[0].set_ylabel('Distance to 2nd Nearest Neighbor (pc)')
axes[1].scatter(all_first, all_third, c='gray', s=3, alpha=0.2, label=r'$\rho$ = '+ str(np.around(scipy.stats.pearsonr(all_first, all_third)[0], 2)))
axes[1].set_xlim((0,1000))
axes[1].set_ylim((0,2000))
axes[1].scatter(xvals2, medians2, c='b')
axes[1].errorbar(xvals2, medians2, c='b', capsize=4, ls='none', xerr=None, yerr=err2)
axes[1].set_title(str(res)+'pc Resolution')
axes[1].set_xlabel('Distance to 1st Nearest Neighbor (pc)')
axes[1].set_ylabel('Distance to 3rd Nearest Neighbor (pc)')
#for i in range(num_bins):
    #axes[0].scatter(xvals[i], medians[i], c='b')
    #axes[0].errorbar(xvals[i], medians[i], c='b', capsize=4, ls='none', xerr=None, yerr=err1)
    #axes[1].scatter(xvals2[i], medians2[i], c='b')
    #axes[1].errorbar(xvals2[i], medians2[i], c='b', capsize=4, ls='none', xerr=None, yerr=[[np.percentile(np.take(all_third, np.where(bins2[2]==(i+1))),16)],[np.percentile(np.take(all_third, np.where(bins2[2]==(i+1))),84)]])

x1 = np.linspace(0,2500, 2501)
axes[0].plot(x1, x1*(np.median(all_second)/np.median(all_first)), c='r', linestyle='--', linewidth=3, label='Slope: Median 2nd / Median 1st')
axes[1].plot(x1, x1*(np.median(all_third)/np.median(all_first)), c='r', linestyle='--', linewidth=3, label='Slope: Median 3rd / Median 1st')
axes[0].legend()
axes[1].legend()
plt.legend()
plt.show()
plt.close()
