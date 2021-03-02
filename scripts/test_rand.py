#######################################################
# Testing random distributions vs. true distributions #
#######################################################


import rand
from astropy.table import Table

image = '/Users/josh/projects/ngc6300_mask.fits'
errImage = '/Users/josh/projects/ngc6300_12m+7m+tp_co21_150pc_strict_emom0.fits'
tab_loc = '/Users/josh/projects/'
t = Table.read(tab_loc+'phangs_sample_table_v1p6.fits')

obj = rand.deprojectMap(image, errImage, t[111]['orient_ra'], t[111]['orient_dec'],
             t[111]['orient_posang'], t[111]['orient_incl'], t[111]['dist'] )

numTests = 5
inten = obj[0]
SNR = obj[2]
dx = obj[5]
dy = obj[6]
galDist = t[111]['dist']
tests = rand.random(numTests, inten, dx, dy, obj[3], obj[4], galDist)

coords = [tests[4], tests[5]]
