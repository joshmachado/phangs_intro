import numpy as np
import scipy.stats
import pandas as pd
import os.path

sources = ['ngc1433', 'ngc2835', 'ngc3621', 'ngc6300']
res = [60, 90, 120, 150]
prop = ['MLUM_MSUN', 'SIGV_KMS']
stats = np.zeros([len(prop), 3])
for i in range(len(res)):
    df = pd.DataFrame(columns = ('source', 'MLUM_1_2', 'MLUM_1_3', 'MLUM_2_3',
    'SIGV_1_2', 'SIGV_1_3', 'SIGV_2_3'))
    for j in range(len(sources)):
        fp = '/Users/josh/projects/intro/'+str(sources[i])+'/'+str(sources[i])+'_'+str(res[j])+'pc_cloud_stats.csv'
        if os.path.isfile(fp) == True:
            cat = pd.read_csv(fp)
            for k in range(len(prop)):
                first = cat[str(prop[k])+'_nn']
                second = cat[str(prop[k])+'_nn2']
                third = cat[str(prop[k])+'_nn3']
                stats[k] = [scipy.stats.spearmanr(first, second)[0], scipy.stats.spearmanr(first, third)[0], scipy.stats.spearmanr(second, third)[0]]
        df.loc[j] = [str(sources[i]), stats[0][0], stats[0][1], stats[0][2], stats[1][0], stats[1][1], stats[1][2]]
    df.to_csv('/Users/josh/projects/intro/stats/spearman_coef_'+str(res[i])+'pc.csv', index=False)
