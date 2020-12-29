### Measuring GMC Property Correlation as Function of Galactocentric Radius###

import numpy as np
import scipy.stats
import pandas as pd
import os.path
import matplotlib.pyplot as plt
from astropy.table import Table



res = [60, 90, 120, 150]
sources = ['ngc1433', 'ngc3621', 'ngc6300', 'ngc2835']
prop = ['MLUM_MSUN', 'SIGV_KMS']#, 'RAD_PC']


for i in range(len(res)):
    for j in range(len(sources)):
        fp = '/Users/josh/projects/intro/'+str(sources[j])+'/'+str(sources[j])+'_'+str(res[i])+'pc_cloud_stats.csv'
        raw = '/Users/josh/projects/intro/'+str(sources[j])+'/'+str(sources[j])+'_12m+7m+tp_co21_'+str(res[i])+'pc_props.fits.bz2'
        if os.path.isfile(fp) == True:
            cat = pd.read_csv(fp)
            tab = Table.read(raw)
            rad_bins = np.linspace(0, np.max(tab['RGAL_KPC']), 5)
            fig, axes = plt.subplots(len(prop), 2, figsize=(10,10))
            for k in range(len(prop)):
                clouds = np.zeros((len(cat), 2))
                clouds_2 = np.zeros((len(cat), 2))
                clouds_3 = np.zeros((len(cat), 2))
                for x in range(len(rad_bins)):
                    for z in range(len(cat)):
                        nn_num = int(cat['nn_index'][z]) -1
                        nn2_num = int(cat['nn2_index'][z]) -1
                        nn3_num = int(cat['nn3_index'][z]) -1

                        if rad_bins[x-1] < tab['RGAL_KPC'][z] < rad_bins[x]:
                            clouds[z] = int(cat['cloudnum'][z])-1, nn_num
                            clouds_2[z] = int(cat['cloudnum'][z])-1, nn2_num
                            clouds_3[z] = int(cat['cloudnum'][z])-1, nn3_num

                    cloud_prop = tab[prop[k]][clouds[clouds[:,0]!=0].astype(int)]
                    prop_nn = tab[prop[k]][clouds[clouds[:,1]!=0][:,1].astype(int)]
                    prop_nn2 = tab[prop[k]][clouds_2[clouds_2[:,1]!=0][:,1].astype(int)]
                    prop_nn3 = tab[prop[k]][clouds_3[clouds_3[:,1]!=0][:,1].astype(int)]

                    #print(str(len(prop_nn)))
                    #print('bin'+str(x))
                    if x!=0:
                        axes[k,0].scatter(prop_nn, prop_nn2[0:len(prop_nn)], alpha=0.35, label='N = '+str(int(len(prop_nn)))+', '+str('%.2f' % rad_bins[x-1])+' < RGAL < '+str('%.2f' % rad_bins[x])+r'kpc, $\rho$ = '+str('%.2f' % scipy.stats.spearmanr(prop_nn, prop_nn2[0:len(prop_nn)], nan_policy='omit')[0]))
                        axes[k,1].scatter(prop_nn2[0:len(prop_nn3)], prop_nn3, alpha=0.35, label='N = '+str(int(len(prop_nn2)))+', '+str('%.2f' % rad_bins[x-1])+' < RGAL < '+str('%.2f' % rad_bins[x])+r'kpc, $\rho$ = '+str('%.2f' % scipy.stats.spearmanr(prop_nn2[0:len(prop_nn)], prop_nn3, nan_policy='omit')[0]))
                        axes[k,0].set_ylabel(prop[k])
                        axes[k,0].legend()
                        axes[k,1].legend()
            axes[0,0].set_title('First to Second')
            axes[0,1].set_title('Second to Third')
            fig.suptitle(sources[j]+' at '+str(res[i])+'pc resolution', fontsize=16)
            plt.show()
            plt.close()
