import scipy.stats
import pandas as pd
from astropy.table import Table
import os.path
import numpy as np

sources = ['ngc1433', 'ngc2835', 'ngc3621', 'ngc6300']
res = [60, 90, 120, 150]
for j in range(len(res)):
    for i in range(len(sources)):
        fp = '/Users/josh/projects/intro/'+str(sources[i])+'/'+str(sources[i])+'_12m+7m+tp_co21_'+str(res[j])+'pc_props.fits.bz2'
        stats = '/Users/josh/projects/intro/'+str(sources[i])+'/'+str(sources[i])+'_'+str(res[j])+'pc_cloud_stats.csv'
        if os.path.isfile(fp) == True:
            tab = Table.read(fp)
            cat = pd.read_csv(stats)
            df = pd.DataFrame(cat)
            prop = ['MLUM_MSUN', 'SIGV_KMS']
            prop_nn = np.zeros(len(cat))
            prop_nn2 = np.zeros(len(cat))
            prop_nn3 = np.zeros(len(cat))
            for z in range(len(prop)):
                for k in range(len(cat)):
                    cloudnum = int(cat['cloudnum'][k])
                    nn = int(cat['nn_index'][cloudnum-1])
                    nn2 = int(cat['nn2_index'][cloudnum-1])
                    nn3 = int(cat['nn3_index'][cloudnum-1])
                    prop_nn[k] = tab[prop[z]][nn-1]
                    prop_nn2[k] = tab[prop[z]][nn2-1]
                    prop_nn3[k] = tab[prop[z]][nn3-1]
                df[str(prop[z])+'_nn'] = prop_nn
                df[str(prop[z])+'_nn2'] = prop_nn2
                df[str(prop[z])+'_nn3'] = prop_nn3
                df.to_csv('/Users/josh/projects/intro/'+str(sources[i])+'/'+str(sources[i])+'_'+str(res[j])+'pc_cloud_stats.csv', index=False)
