from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import pandas as pd
import os.path

sources = ['ngc6300', 'ngc3621']
loc = '/Users/josh/projects/intro/'
res = np.array([60,90,120,150])
for i in range(len(res)):
    for j in range(len(sources)):
        fp = '/Users/josh/projects/intro/'+sources[j]+'/max/'+sources[j]+'_12m+7m+tp_co21_'+str(res[i])+'pc_props.fits.bz2'
        if os.path.isfile(fp) == True:
            tab = Table.read(fp)

            x_ctr = np.array(tab['XCTR_DEG'])
            x_ctr = x_ctr*np.cos(np.deg2rad(x_ctr))
            y_ctr = np.array(tab['YCTR_DEG'])

            x_max = np.array(tab['XMAX_DEG'])
            x_max = x_max*np.cos(np.deg2rad(x_max))
            y_max = np.array(tab['YMAX_DEG'])

            plt.scatter(x_ctr, y_ctr, c='red', alpha=0.5, label='Centers')
            plt.scatter(x_max, y_max, c='blue', alpha=0.5, label='Max')
            plt.xlabel('R.A.')
            plt.ylabel('DEC')
            plt.legend()
            plt.title(str(sources[j])+' ('+str(res[i])+'pc_homogenized) - Maxes vs. Centers', fontsize=15)
            plt.savefig(str(sources[j])+'_max_v_center_loc'+str(res[i])+'pc.pdf')
            plt.close()
        else:
            print(str(sources[j])+' '+str(res[i])+'pc not found')
