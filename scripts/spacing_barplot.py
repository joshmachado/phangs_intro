from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import pandas as pd
import matplotlib.lines as mlines
import os.path


loc = '/Users/josh/projects/intro/sources/'
res = 150
sources = pd.read_csv('sources.csv')
sources = sources['circinus'].values.tolist()
center_first, center_second, center_third = np.zeros(len(sources)), np.zeros(len(sources)), np.zeros(len(sources))
upper_first, upper_second, upper_third = np.zeros(len(sources)), np.zeros(len(sources)), np.zeros(len(sources))
lower_first, lower_second, lower_third = np.zeros(len(sources)), np.zeros(len(sources)), np.zeros(len(sources))
all_first, all_second, all_third = [], [], []
for i in range(len(sources)):
    fp = loc+str(sources[i])+'/'+str(sources[i])+'_'+str(res)+'pc_cloud_stats.csv'
    if os.path.isfile(fp):
        cat = pd.read_csv(fp)

        #EACH INDIVIDUAL GALAXY
        center_first[i] = np.percentile(cat['min_dist'].tolist(), 50)
        upper_first[i] = np.percentile(cat['min_dist'].tolist(), 84)
        lower_first[i] = np.percentile(cat['min_dist'].tolist(), 16)

        center_second[i] = np.percentile(cat['min_dist2nd'].tolist(), 50)
        upper_second[i] = np.percentile(cat['min_dist2nd'].tolist(), 84)
        lower_second[i] = np.percentile(cat['min_dist2nd'].tolist(), 16)

        center_third[i] = np.percentile(cat['min_dist3rd'].tolist(), 50)
        upper_third[i] = np.percentile(cat['min_dist3rd'].tolist(), 84)
        lower_third[i] = np.percentile(cat['min_dist3rd'].tolist(), 16)

        #ALL GALAXIES
        all_first = all_first + cat['min_dist'].tolist()
        all_second = all_second + cat['min_dist2nd'].tolist()
        all_third = all_third + cat['min_dist3rd'].tolist()


center_first = center_first[center_first != 0]
center_second = center_second[center_second != 0]
center_third = center_third[center_third != 0]
lower_first = lower_first[lower_first != 0]
lower_second = lower_second[lower_second != 0]
lower_third = lower_third[lower_third != 0]
upper_first = upper_first[upper_first != 0]
upper_second = upper_second[upper_second != 0]
upper_third = upper_third[upper_third != 0]
err1 = np.abs([lower_first, upper_first] - center_first)
err2 = np.abs([lower_second, upper_second] - center_second)
err3 = np.abs([lower_third, upper_third] - center_third)

if len(center_first) % 2 == 1:
    xs = np.concatenate((np.linspace(0, int(len(center_first) / 2), int(len(center_first) / 2)+1),
    np.linspace(int(len(center_first) / 2)+6, len(center_first)+4, int(len(center_first) / 2))))

else:
    xs = np.concatenate((np.linspace(0,len(center_first)/2,int(len(center_first)/2) +1),np.linspace(len(center_first)/2 +6,len(center_first)+4,int(len(center_first)/2)-1)))

alpha = 0.6
plt.scatter(xs, center_first, label='1st Nearest Neighbor Separation (pc)', alpha = alpha)
plt.errorbar(xs, center_first, xerr=None, yerr=err1, capsize=4, ls='none', alpha = alpha)

plt.scatter(xs, center_second, label='2nd Nearest Neighbor Separation (pc)', alpha = alpha)
plt.errorbar(xs, center_second, xerr=None, yerr=err2, capsize=4, ls='none', alpha = alpha)

plt.scatter(xs, center_third, label='3rd Nearest Neighbor Separation (pc)', alpha = alpha)
plt.errorbar(xs, center_third, xerr=None, yerr=err3, capsize=4, ls='none', alpha = alpha)

#GLOBAL Median
x = int(len(center_first)/2)+3
plt.scatter(x, np.percentile(all_first, 50), label='All Sources - 1st Nearest', c='red')
plt.errorbar(x, np.percentile(all_first, 50), xerr=None, yerr=[[np.percentile(all_first, 16)], [np.percentile(all_first, 84)]], capsize=4, ls='none', c='red')

plt.scatter(x, np.percentile(all_second, 50), label='All Sources - 2nd Nearest', c='black')
plt.errorbar(x, np.percentile(all_second, 50), xerr=None, yerr=[[np.percentile(all_second, 16)], [np.percentile(all_second, 84)]], capsize=4, ls='none', c='black')

plt.scatter(x, np.percentile(all_third, 50), label='All Sources - 3rd Nearest', c='gray')
plt.errorbar(x, np.percentile(all_third, 50), xerr=None, yerr=[[np.percentile(all_third, 16)], [np.percentile(all_third, 84)]], capsize=4, ls='none', c='gray')

plt.legend()
plt.ylim((0,1500))
plt.ylabel('Separation (pc)')
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False)
plt.title(str(res)+' pc Resolution', fontsize=16)
plt.show()
plt.close()
