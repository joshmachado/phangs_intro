import numpy as np
import scipy.stats
import pandas as pd
import os.path

sources = ['ngc1433', 'ngc2835', 'ngc3621', 'ngc6300']
res = [60, 90, 120, 150]
for j in range(len(res)):
    df = pd.DataFrame(columns = ('source', '1_2', '1_3', '2_3'))
    for i in range(len(sources)):
        fp = '/Users/josh/projects/intro/'+str(sources[i])+'/'+str(sources[i])+'_'+str(res[j])+'pc_cloud_stats.csv'
        if os.path.isfile(fp) == True:
            cat = pd.read_csv(fp)
            first = cat['beam_sep_nn']
            second = cat['beam_sep_nn2']
            third = cat['beam_sep_nn3']
            df.loc[i] = [str(sources[i]), scipy.stats.pearsonr(first, second)[0], scipy.stats.pearsonr(first, third)[0], scipy.stats.pearsonr(second, third)[0]]
    df.to_csv('/Users/josh/projects/intro/stats/pearson_coef_'+str(res[j])+'pc.csv', index=False)
