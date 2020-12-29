from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import pandas as pd

fp = '/Users/josh/projects/intro/ngc3621_12m+7m+tp_co21_90pc_props.fits.bz2'
#fp = '/data/kant/0/sun.1608/PHANGS/ALMA/v3p4-CPROPS/STv1p5/90pc_homogenized/ngc3621_12m+7m+tp_co21_90pc_props.fits.bz2'
cat = pd.read_csv('ngc3621_sources.csv')

tab = Table.read(fp)
x = np.array(tab['XCTR_DEG'])
y = np.array(tab['YCTR_DEG'])
ys = y*np.cos(x)
vel = tab['VCTR_KMS']

plt.scatter(x,y, c=vel, cmap='RdBu', alpha=0.75)
plt.xlabel('R.A.')
plt.ylabel('DEC')
plt.title('NGC3621 (90pc_homogenized)', fontsize=15)
plt.colorbar(label='Velocity km/s')
plt.savefig('ngc3621_peaks.pdf')
plt.close()

plt.scatter(x,ys, c=vel, cmap='RdBu', alpha=0.75)
plt.xlabel('R.A.')
plt.ylabel('DEC')
plt.title('NGC3621 (90pc_homogenized) - Corrected', fontsize=15)
plt.colorbar(label='Velocity km/s')
plt.savefig('ngc3621_peaks_corr.pdf')
plt.close()

plt.hist(cat.true_dist, bins = (np.linspace(2000, 10000, 10)), label=r'Eucledian (N = '+str(len(x))+')', alpha = 0.5)
plt.hist(cat.true_dist_corr, bins = (np.linspace(2000, 10000, 10)), label=r'Spherical Corrected (N ='+str(len(x))+')', alpha = 0.5)
plt.xlabel('Distance to Nearest Neighbor (pc)')
plt.ylabel('Counts')
plt.legend()
plt.title('NGC3621 (90pc_homogenized)')
plt.savefig('ngc3621_hist.pdf')
plt.close()

#fig, axes = plt.subplots(1,2, figsize=(17,7))

#fig1 = axes[0].scatter(x, y, c=vel, cmap='RdBu', alpha = 0.5, label='Eucledian')
#fig2 = axes[1].scatter(x, ys, c=vel, cmap='RdBu',alpha = 0.5, label='Angle Corrected')
#axes[0].set_xlabel('R.A.')
#axes[0].set_ylabel('DEC')
#axes[0].set_label('HEY')
#axes[1].set_xlabel('R.A.')
#axes[1].set_ylabel('DEC')
#cbar1 = fig.colorbar(fig1,ax=axes[0], cmap='RdBu')
#cbar2 = fig.colorbar(fig1,ax=axes[1], cmap='RdBu')
#cbar1.set_label('Velocity km/s')
#cbar2.set_label('Velocity km/s')

#fig2 = axes[1,0].scatter(x, y, c=vel, cmap='RdBu', alpha = 0.5, label='Eucledian')
#axes[1,0].set_xlabel('R.A.')
#axes[1,0].set_ylabel('DEC')
#cbar = fig.colorbar(fig1,ax=axes[1,0], cmap='RdBu')
#cbar.set_label('Velocity km/s')

#fig3 = axes[1,1].scatter(x, ys, c=vel, cmap='RdBu',alpha = 0.5, label='Angle Corrected')
#axes[1,1].set_xlabel('R.A.')
#axes[1,1].set_ylabel('DEC')
#cbar = fig.colorbar(fig1,ax=axes[1,1], cmap='RdBu')
#cbar.set_label('Velocity km/s')


#fig3 = plt.hist(cat.true_dist, bins = (np.linspace(2000, 10000, 10)), label='Eucledian', alpha = 0.5)
#fig3.hist(cat.true_dist_corr, bins = (np.linspace(2000, 10000, 10)), label='Distance Corrected', alpha = 0.5)
#axes[0,1].set_xlabel('Distance to Nearest Neighbor (pc)')
#axes[0,1].set_ylabel('Counts')
#fig.suptitle('NGC3621 (90pc_homogenized)', fontsize=18)
#plt.show()
