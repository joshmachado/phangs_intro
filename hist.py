from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import pandas as pd

#fp = '/Users/josh/projects/intro/ngc3621_12m+7m+tp_co21_90pc_props.fits.bz2'
#fp = '/data/kant/0/sun.1608/PHANGS/ALMA/v3p4-CPROPS/STv1p5/90pc_homogenized/ngc3621_12m+7m+tp_co21_90pc_props.fits.bz2'
res = np.array([60,90,120,150])

cat60 = pd.read_csv('ngc3621_60pc.csv')
cat90 = pd.read_csv('ngc3621_90pc.csv')
cat120 = pd.read_csv('ngc3621_120pc.csv')
cat150 = pd.read_csv('ngc3621_150pc.csv')

plt.hist(cat60.true_dist, bins = (np.linspace(2000, 10000, 10)), label=r'60pc (N = '+str(len(cat60.true_dist))+')', alpha = 0.25)
plt.hist(cat90.true_dist, bins = (np.linspace(2000, 10000, 10)), label=r'90pc (N = '+str(len(cat90.true_dist))+')', alpha = 0.25)
plt.hist(cat120.true_dist, bins = (np.linspace(2000, 10000, 10)), label=r'120pc (N = '+str(len(cat120.true_dist))+')', alpha = 0.25)
plt.hist(cat150.true_dist, bins = (np.linspace(2000, 10000, 10)), label=r'150pc (N = '+str(len(cat150.true_dist))+')', alpha = 0.25)
plt.legend()
plt.xlabel('Distance to Nearest Neighbor (pc)')
plt.ylabel('Counts')
plt.title('NGC3621', fontsize=15)
plt.savefig('hist_all.pdf')
