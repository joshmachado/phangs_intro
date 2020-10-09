from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import pandas as pd
loc = '/Users/josh/projects/intro'
source = 'ngc3621'
res = np.array([60,90,120,150])
for i in range(len(res)):

    tab = Table.read(str(source)+'_12m+7m+tp_co21_'+str(res[i])+'pc_props.fits.bz2')
    cat = pd.read_csv(str(source)+'_'+str(res[i])+'pc_cloud_stats.csv')
    inner = np.zeros(len(cat))
    outer = np.zeros(len(cat))
    cutoff = 3.5
    #cutoff = np.median(tab['RGAL_KPC'])
    for j in range(len(tab)):
        if tab['RGAL_KPC'][j] < cutoff:
            inner[j] = cat.min_dist3rd[j]
        else:
            outer[j] = cat.min_dist3rd[j]

    inner = inner[inner != 0]
    outer = outer[outer != 0]
    plt.hist(inner, bins = (np.linspace(res[i], 600, 15)), label='Inner', alpha=0.5)
    plt.hist(outer, bins = (np.linspace(res[i], 600, 15)), label='Outer', alpha=0.5)
    plt.legend()
    plt.title('NGC3621 at '+str(res[i])+'pc resolution - Cutoff = '+str(np.round(cutoff, 2))+'kpc')
    plt.savefig(loc+'/plots/inner_outer_'+str(res[i])+'pc')
    plt.close()
