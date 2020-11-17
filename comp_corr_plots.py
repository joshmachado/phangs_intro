import numpy as np
import scipy.stats
import pandas as pd
import os.path
import matplotlib.pyplot as plt

sources = ['ngc1433', 'ngc3621', 'ngc6300', 'ngc2835']
res = [60, 90, 120, 150]
prop = ['MLUM_MSUN', 'SIGV_KMS', 'RAD_PC']



for i in range(len(res)):
    for j in range(len(sources)):
        fp = '/Users/josh/projects/intro/'+str(sources[j])+'/'+str(sources[j])+'_'+str(res[i])+'pc_cloud_stats.csv'
        if os.path.isfile(fp) == True:
            cat = pd.read_csv(fp)
            fig, axes = plt.subplots(3, 2, figsize=(10,10))
            for k in range(len(prop)):
                first = np.array(cat[str(prop[k])+'_nn'])
                second = np.array(cat[str(prop[k])+'_nn2'])
                third = np.array(cat[str(prop[k])+'_nn3'])
                axes[k,0].scatter(first, second)
                axes[k,1].scatter(second, third)
                axes[k,0].set_ylabel(prop[k])
            axes[0,0].set_title('First to Second')
            axes[0,1].set_title('Second to Third')
            fig.suptitle(sources[j]+' at '+str(res[i])+'pc resolution', fontsize=16)
            plt.savefig('/Users/josh/projects/intro/stats/'+sources[j]+'_'+str(res[i])+'_propcorr.pdf')
