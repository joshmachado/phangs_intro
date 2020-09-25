from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import pandas as pd
import aplpy
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import numpy as np
from astropy.coordinates import SkyCoord
loc = '/Users/josh/projects/intro'

m0 = fits.open(loc+'/maps/ngc3621_1_12m+7m_co21_60pc_broad_mom0.fits')[0]
sources = pd.read_csv('ngc3621_60pc.csv')


mom0  = aplpy.FITSFigure(m0.data)
mom0.show_colorscale()
mom0.add_colorbar()
mom0.colorbar.set_axis_label_font(size=16)
mom0.colorbar.set_axis_label_text('K')
#tkin.colorbar.set_font(size=14)
#tkin.add_scalebar(0.015)
#tkin.scalebar.set_label('80pc')
#mom0.add_beam()
#tkin.set_title('Kinematic Temperature')
#tkin.axis_labels.set_font(size=19)
#tkin.tick_labels.set_font(size=14)





base_wcs = WCS(m0.header)
ax = plt.subplot(projection=base_wcs.celestial)
ax.imshow(m0.data, origin='lower')
coords = SkyCoord(ra = sources.corr_x*u.degree, dec=sources.y*u.degree)
ax.scatter(coords.ra, coords.dec, facecolors='none', edgecolors='r')
plt.show()
