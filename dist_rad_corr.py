import pandas as pd
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats

rad_bins = [0,1,2,3,4,5]
cat = pd.read_csv('/Users/josh/projects/intro/ngc3621/ngc3621_60pc_cloud_stats.csv')
cloudnum = np.array(cat['cloudnum'])
first = np.array(cat['min_dist'])
second = np.array(cat['min_dist2nd'])
third = np.array(cat['min_dist3rd'])

tab = Table.read('/Users/josh/projects/intro/ngc3621/ngc3621_12m+7m+tp_co21_60pc_props.fits.bz2')
rad = np.array(tab['RGAL_KPC']).byteswap().newbyteorder()


df = pd.DataFrame({'cloudnum':cloudnum, 'first':first, 'second':second,
'third':third, 'RGAL_KPC':rad})

df = df.sort_values(by=['RGAL_KPC'])
radregion = []

for i in range(len(rad_bins)):
    radregion.append(df[df['RGAL_KPC'].between(rad_bins[i-1],rad_bins[i], inclusive=False)])
radregion = radregion[1:len(radregion)]


fig_size = (7.5,8.5)
fig1, ax1 = plt.subplots(len(radregion), 1, figsize=fig_size, dpi=175) #Separation v galactocentric radius
if i!=0:
    for k in range(len(radregion)):

        ax1[k].scatter(radregion[k]['first'], radregion[k]['RGAL_KPC'], color='b', marker='.', alpha=0.2, label='N = '+str(int(len(radregion[k]['first'])))+', '+str('%.2f' % rad_bins[k])+' < RGAL < '+str('%.2f' % rad_bins[k+1])+r'kpc, $\rho$ = '+str('%.2f' % scipy.stats.spearmanr(radregion[k]['first'], radregion[k]['RGAL_KPC'], nan_policy='omit')[0]))
        ax1[k].scatter(radregion[k]['second'], radregion[k]['RGAL_KPC'], color='r', marker='v', alpha=0.2, label='N = '+str(int(len(radregion[k]['second'])))+', '+str('%.2f' % rad_bins[k])+' < RGAL < '+str('%.2f' % rad_bins[k+1])+r'kpc, $\rho$ = '+str('%.2f' % scipy.stats.spearmanr(radregion[k]['second'], radregion[k]['RGAL_KPC'], nan_policy='omit')[0]))
        ax1[k].scatter(radregion[k]['third'], radregion[k]['RGAL_KPC'], color='g', marker='s', alpha=0.2, label='N = '+str(int(len(radregion[k]['third'])))+', '+str('%.2f' % rad_bins[k])+' < RGAL < '+str('%.2f' % rad_bins[k+1])+r'kpc, $\rho$ = '+str('%.2f' % scipy.stats.spearmanr(radregion[k]['third'], radregion[k]['RGAL_KPC'], nan_policy='omit')[0]))
        ax1[k].legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
        ax1[k].get_xaxis().set_visible(False)


ax1[0].set_title('Separation vs. Galactocentric Radius')
ax1[k].set_xlabel('Separation')
ax1[k].get_xaxis().set_visible(True)
plt.tight_layout()
fig1.text(0.01,0.5, 'RGAL_KPC', va='center', rotation='vertical')
#fig1.suptitle('NGC 3621', fontsize=16)

plt.show()
plt.close()
