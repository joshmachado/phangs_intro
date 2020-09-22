from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import pandas as pd

#fp = '/Users/josh/projects/intro/ngc3621_12m+7m+tp_co21_90pc_props.fits.bz2'
#fp = '/data/kant/0/sun.1608/PHANGS/ALMA/v3p4-CPROPS/STv1p5/90pc_homogenized/ngc3621_12m+7m+tp_co21_90pc_props.fits.bz2'
res = np.array([60,90,120,150])
for i in range(len(res)):
    cat = pd.read_csv('ngc3621_'+str(res[i])+'pc.csv')
    fp = '/Users/josh/projects/intro/ngc3621_12m+7m+tp_co21_'+str(res[i])+'pc_props.fits.bz2'
    tab = Table.read(fp)
    x = np.array(tab['XCTR_DEG'])
    y = np.array(tab['YCTR_DEG'])
    ys = y*np.cos(x)
    vel = tab['VCTR_KMS']

    plt.scatter(x,y, c=vel, cmap='RdBu', alpha=0.75)
    plt.xlabel('R.A.')
    plt.ylabel('DEC')
    plt.title('NGC3621 ('+str(res[i])+'pc_homogenized)', fontsize=15)
    plt.colorbar(label='Velocity km/s')
    plt.savefig('ngc3621_peaks'+str(res[i])+'pc.pdf')
    plt.close()

    plt.scatter(x,ys, c=vel, cmap='RdBu', alpha=0.75)
    plt.xlabel('R.A.')
    plt.ylabel('DEC')
    plt.title('NGC3621 ('+str(res[i])+'pc_homogenized) - Corrected', fontsize=15)
    plt.colorbar(label='Velocity km/s')
    plt.savefig('ngc3621_peaks_corr'+str(res[i])+'pc.pdf')
    plt.close()

    plt.hist(cat.true_dist, bins = (np.linspace(2000, 10000, 10)), label=r'Eucledian (N = '+str(len(x))+')', alpha = 0.5)
    plt.hist(cat.true_dist_corr, bins = (np.linspace(2000, 10000, 10)), label=r'Spherical Corrected (N ='+str(len(x))+')', alpha = 0.5)
    plt.xlabel('Distance to Nearest Neighbor (pc)')
    plt.ylabel('Counts')
    plt.legend()
    plt.title('NGC3621 ('+str(res[i])+'pc_homogenized)')
    plt.savefig('ngc3621_hist'+str(res[i])+'pc.pdf')
    plt.close()
