import numpy as np
import scipy.stats
import pandas as pd
import os.path
import matplotlib.pyplot as plt

sources = ['ngc1433', 'ngc3621', 'ngc6300', 'ngc2835']
res = [60, 90, 120, 150]
prop = ['MLUM_MSUN', 'SIGV_KMS', 'RAD_NODC']
stats = np.zeros([len(prop), 6])
for i in range(len(res)):
    df = pd.DataFrame(columns = ('source', 'first_second', 'pval_first_second', 'first_third', 'pval_first_third',
    'second_third', 'pval_second_third'))
    for j in range(len(sources)):
        fp = '/Users/josh/projects/intro/'+str(sources[j])+'/'+str(sources[j])+'_'+str(res[i])+'pc_cloud_stats.csv'
        if os.path.isfile(fp) == True:
            cat = pd.read_csv(fp)
            for k in range(len(prop)):
                first = np.array(cat[str(prop[k])+'_nn'])
                second = np.array(cat[str(prop[k])+'_nn2'])
                third = np.array(cat[str(prop[k])+'_nn3'])
                fs = scipy.stats.spearmanr(first, second, nan_policy='omit')
                ft = scipy.stats.spearmanr(first, third, nan_policy='omit')
                st = scipy.stats.spearmanr(second, third, nan_policy='omit')
                stats[k] = [fs[0], fs[1], ft[0], ft[1], st[0], st[1]]
                df.loc[j] = [str(sources[j]), stats[0][0], stats[0][1], stats[0][2], stats[0][3], stats[0][4], stats[0][5]]
                df.to_csv('/Users/josh/projects/intro/stats/'+str(prop[k])+'_spearman_coef_'+str(res[i])+'pc.csv', index=False)
                print('Finished '+str(res[i])+' '+sources[j]+' '+prop[k])
