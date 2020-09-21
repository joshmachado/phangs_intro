from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm

fp = '/Users/josh/projects/intro/ngc3621_12m+7m+tp_co21_90pc_props.fits.bz2'
#fp = '/data/kant/0/sun.1608/PHANGS/ALMA/v3p4-CPROPS/STv1p5/90pc_homogenized/ngc3621_12m+7m+tp_co21_90pc_props.fits.bz2'


tab = Table.read(fp)
x = np.array(tab['XCTR_DEG'])
y = np.array(tab['YCTR_DEG'])
ys = y*np.cos(x)
vel = tab['VCTR_KMS']
plt.scatter(x, y, c=vel, cmap='RdBu')
plt.xlabel('R.A.')
plt.ylabel('DEC')
cbar = plt.colorbar(cmap='RdBu')
cbar.set_label('Velocity km/s')

#plt.colorbar()
plt.show()
