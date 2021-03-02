from astropy.table import Table
import numpy as np
import peak_stats
import pandas as pd
import matplotlib.pyplot as plt
import os.path
import matplotlib.patches as mpatches
import scipy.stats


loc = '/Users/josh/projects/intro/sources/'
res = 150
res_pc = [60,90,120,150]
sources = pd.read_csv('sources.csv')
sources = sources['circinus'].values.tolist()

tab_loc = '/Users/josh/projects/'
t = Table.read(tab_loc+'phangs_sample_table_v1p6.fits')
names = list(t['name'])
#names = set(sources) & set(names)
lco = t['lco_phangs']
sfr = t['props_sfr']
m_star = t['props_mstar']
size_reff = t['size_reff']
median, sfr_val, m_star_val, size_reff_val = np.zeros(len(sources)), np.zeros(len(sources)), np.zeros(len(sources)), np.zeros(len(sources))
lco_val = np.zeros(len(sources))
dist = np.zeros(len(sources))
err_bar = np.zeros((2, len(sources)))
for i in range(len(sources)):
    fp = loc+sources[i]+'/'+sources[i]+'_stats.csv'
    if os.path.isfile(fp):
        stats = pd.read_csv(fp)
        median[i] = stats['50th_nn'][res_pc.index(res)]
        err_bar[0][i] = stats['16th_nn'][res_pc.index(res)]
        err_bar[1][i] = stats['84th_nn'][res_pc.index(res)]
        sfr_val[i] = sfr[names.index(sources[i])]
        lco_val[i] = lco[names.index(sources[i])]
        m_star_val[i] = m_star[names.index(sources[i])]
        size_reff_val[i] = size_reff[names.index(sources[i])]
        dist[i] = t['dist'][names.index(sources[i])]

#Remove zero median Sources
ind = np.where(median != 0)
median = median[ind]
lco_val = lco_val[ind]
sfr_val = sfr_val[ind]
m_star_val = m_star_val[ind]
size_reff_val = size_reff_val[ind]
dist = dist[ind]
err = np.zeros((2, len(median)))
err[0] = err_bar[0][err_bar[0] != 0]
err[1] = err_bar[1][err_bar[1] != 0]

#arcsec to parsec
size_reff_val = size_reff_val *np.pi/(180 * 3600)* dist * 1e3


#props per area (per pc^2)
lco_area = lco_val / (np.pi * (size_reff_val*1000)**2)
sfr_area = sfr_val / (np.pi * (size_reff_val*1000)**2)
m_star_area = m_star_val / (np.pi * (size_reff_val*1000)**2)
##plots

figsize = (15,10)
fig, ax = plt.subplots(2, 2, figsize=figsize)

ax[0,0].scatter(sfr_val, median, edgecolor='black', label=r'Spearman: $\rho$ = '+ str(np.around(scipy.stats.spearmanr(sfr_val, median)[0], 2)))
ax[0,0].errorbar(sfr_val, median, xerr=None, yerr=err, capsize=4, ls='none', c='gray', alpha=0.3)
ax[0,0].set_title('SFR v. Median Nearest Neighbor Separation')
ax[0,0].set_xlabel('SFR [M'+r'$_{\odot}$/yr]')
ax[0,0].set_ylabel('Separation [pc]')
ax[0,0].set_xscale('log')
ax[0,0].set_ylim((0,650))
ax[0,0].legend()

ax[0,1].scatter(m_star_val, median, edgecolor='black', label=r'Spearman: $\rho$ = '+ str(np.around(scipy.stats.spearmanr(m_star_val, median)[0], 2)))
ax[0,1].errorbar(m_star_val, median, xerr=None, yerr=err, capsize=4, ls='none', c='gray', alpha=0.3)
ax[0,1].set_xscale('log')
ax[0,1].set_title('Mass v. Median Nearest Neighbor Separation')
ax[0,1].set_xlabel('Mass [M'+r'$_{\odot}$]')
ax[0,1].set_ylim((0,650))
ax[0,1].legend()

ax[1,0].scatter(size_reff_val, median, edgecolor='black', label=r'Spearman: $\rho$ = '+ str(np.around(scipy.stats.spearmanr(size_reff_val, median)[0], 2)))
ax[1,0].errorbar(size_reff_val, median, xerr=None, yerr=err, capsize=4, ls='none', c='gray', alpha=0.3)
ax[1,0].set_title(r'R$_{eff}$'+' v. Median Nearest Neighbor Separation')
ax[1,0].set_xlabel(r'R$_{eff}$'+' [kpc]')
ax[1,0].set_ylim((0,650))
ax[1,0].legend()

ax[1,1].scatter(lco_val, median, edgecolor='black', label=r'Spearman: $\rho$ = '+ str(np.around(scipy.stats.spearmanr(lco_val, median)[0], 2)))
ax[1,1].errorbar(lco_val, median, xerr=None, yerr=err, capsize=4, ls='none', c='gray', alpha=0.3)
ax[1,1].set_title(r'L$_{CO}$'+' v. Median Nearest Neighbor Separation')
ax[1,1].set_xlabel(r'L$_{CO}$'+r' [K km pc$^2$ / s]')
ax[1,1].set_xscale('log')
ax[1,1].set_ylim((0,650))
ax[1,1].legend()

plt.show()
plt.close()

figsize = (15,5)
fig, ax = plt.subplots(1, 3, figsize=figsize)

ax[0].scatter(sfr_area, median, edgecolor='black', label=r'Spearman: $\rho$ = '+ str(np.around(scipy.stats.spearmanr(sfr_area, median)[0], 2)))
ax[0].errorbar(sfr_area, median, xerr=None, yerr=err, capsize=4, ls='none', c='gray', alpha=0.3)
ax[0].set_title('SFR per area v. Median Nearest Neighbor Separation')
ax[0].set_xlabel('SFR per area [M'+r'$_{\odot}$/yr / pc$^2$]')
ax[0].set_ylabel('Separation [pc]')
ax[0].set_xscale('log')
ax[0].set_ylim((0,650))
ax[0].legend()

ax[1].scatter(m_star_area, median, edgecolor='black', label=r'Spearman: $\rho$ = '+ str(np.around(scipy.stats.spearmanr(m_star_area, median)[0], 2)))
ax[1].errorbar(m_star_area, median, xerr=None, yerr=err, capsize=4, ls='none', c='gray', alpha=0.3)
ax[1].set_xscale('log')
ax[1].set_title('Mass per area v. Median Nearest Neighbor Separation')
ax[1].set_xlabel('Mass per area [M'+r'$_{\odot}$ / pc$^2$]')
ax[1].set_ylim((0,650))
ax[1].legend()

ax[2].scatter(lco_area, median, edgecolor='black', label=r'Spearman: $\rho$ = '+ str(np.around(scipy.stats.spearmanr(lco_area, median)[0], 2)))
ax[2].errorbar(lco_area, median, xerr=None, yerr=err, capsize=4, ls='none', c='gray', alpha=0.3)
ax[2].set_title(r'L$_{CO}$ per area'+' v. Median Nearest Neighbor Separation')
ax[2].set_xlabel(r'L$_{CO}$ per area'+r' [K km / s]')
ax[2].set_xscale('log')
ax[2].set_ylim((0,650))
ax[2].legend()
plt.tight_layout()
plt.show()
