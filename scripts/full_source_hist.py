import peak_stats
import pandas as pd
import matplotlib.pyplot as plt
import os.path
import numpy as np
import matplotlib.patches as mpatches


loc = '/Users/josh/projects/intro/sources/'
res = 150
sources = pd.read_csv('sources.csv')
sources = sources['circinus'].values.tolist()
first_sep = []
second_sep = []
third_sep = []
for i in range(len(sources)):
    fp = loc+str(sources[i])+'/'+str(sources[i])+'_'+str(res)+'pc_cloud_stats.csv'
    if os.path.isfile(fp):
        cat = pd.read_csv(fp)
        first_sep = first_sep + cat['min_dist'].tolist()
        second_sep = second_sep + cat['min_dist2nd'].tolist()
        third_sep = third_sep + cat['min_dist3rd'].tolist()

bins =  np.linspace(0,500,50)
log = False
plt.hist(first_sep, bins=bins, lw=1, fc=(0, 0, 0, 0.35), log=log,
label='First, Median Separation: '+str(np.around(np.median(first_sep), decimals=0))+' pc', hatch='/', histtype='step', fill=True, edgecolor='black')
#plt.title('Distance to First Nearest Neighbor (pc)', fontsize=16)
plt.xlabel('pc')
plt.legend()
#plt.savefig('/Users/josh/projects/intro/plots/full_sep_first.pdf')
#plt.close()

plt.hist(second_sep, bins=bins, log=log, ls='dotted', lw=3, fc=(1, 0, 0, 0.5), histtype='step',
label='Second, Median Separation: '+str(np.around(np.median(second_sep), decimals=0))+' pc')
#plt.title('Distance to Second Nearest Neighbor (pc)', fontsize=16)
plt.xlabel('pc')
plt.legend()
#plt.savefig('/Users/josh/projects/intro/plots/full_sep_second.pdf')
#plt.close()

plt.hist(third_sep, bins=bins, log=log, ls='dashed', lw=3, fc=(0, 0, 1, 0.5), histtype='step',
label='Third, Median Separation: '+str(np.around(np.median(third_sep), decimals=0))+' pc')
plt.title(str(res)+'pc Res: Distance to Nearest Neighbors (pc)', fontsize=16)
plt.xlabel('pc')
plt.legend()
#plt.savefig('/Users/josh/projects/intro/plots/full_sep_third.pdf')
#plt.close()
#plt.savefig('/Users/josh/projects/intro/plots/zoomed_linear_'+str(res)+'.pdf')
plt.show()
plt.close()
