from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import pandas as pd
loc = '/Users/josh/projects/intro'

#fp = '/Users/josh/projects/intro/ngc3621_12m+7m+tp_co21_90pc_props.fits.bz2'
#fp = '/data/kant/0/sun.1608/PHANGS/ALMA/v3p4-CPROPS/STv1p5/90pc_homogenized/ngc3621_12m+7m+tp_co21_90pc_props.fits.bz2'
res = np.array([60,90,120,150])
for i in range(len(res)):
    cat = pd.read_csv('ngc3621_'+str(res[i])+'pc.csv')
    fp = '/Users/josh/projects/intro/ngc3621_12m+7m+tp_co21_'+str(res[i])+'pc_props.fits.bz2'
    tab = Table.read(fp)
    x = np.array(tab['XCTR_DEG'])
    y = np.array(tab['YCTR_DEG'])
    xs = x*np.cos(np.deg2rad(x))
    vel = tab['VCTR_KMS']

    #plt.scatter(x,y, c=vel, cmap='RdBu', alpha=0.75)
    #plt.xlabel('DEC')
    #plt.ylabel('R.A.')
    #plt.title('NGC3621 ('+str(res[i])+'pc_homogenized)', fontsize=15)
    #plt.colorbar(label='Velocity km/s')
    #plt.savefig('ngc3621_peaks'+str(res[i])+'pc.pdf')
    #plt.close()

    plt.scatter(xs,y, c=vel, cmap='RdBu', alpha=0.75)
    plt.xlabel('R.A.')
    plt.ylabel('DEC')
    plt.title('NGC3621 ('+str(res[i])+'pc_homogenized) - Dist. Corrected', fontsize=15)
    plt.colorbar(label='Velocity km/s')
    plt.savefig(loc+'/plots/ngc3621_peaks'+str(res[i])+'pc.pdf')
    plt.close()
    print(str(res[i])+'pc peak')

    plt.hist(cat.true_dist_corr, bins = (np.linspace(0, 300, 10)), label=r'Eucledian (N = '+str(len(x))+')', alpha = 0.5)
    #14 sources above 500pc
    plt.xlabel('Distance to Nearest Neighbor (pc)')
    plt.ylabel('Counts')
    plt.legend()
    plt.title('NGC3621 ('+str(res[i])+'pc_homogenized)')
    plt.savefig(loc+'/plots/ngc3621_hist'+str(res[i])+'pc.pdf')
    plt.close()
    print(str(res[i])+'pc hist')
