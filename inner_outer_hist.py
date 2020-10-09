from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import pandas as pd
loc = '/Users/josh/projects/intro'




want to plot separations based on inner / outer

tab = Table.read('ngc3621_12m+7m+tp_co21_150pc_props.fits.bz2')
cat = pd.read_csv('ngc3621_150pc_cloud_stats.csv')
inner = np.zeros(len(cat))
outer = np.zeros(len(cat))
for i in range(len(x)):
    if tab['RGAL_KPC'][i] < np.median(tab['RGAL_KPC']):
        inner[i] = cat.min_dist3rd[i]
    else:
        outer[i] = cat.min_dist3rd[i]

inner = inner[inner != 0]
outer = outer[outer != 0]
plt.hist(inner, bins = (np.linspace(0, 400, 15)), label='Inner', alpha=0.5)
plt.hist(outer, bins = (np.linspace(0, 400, 15)), label='Outer', alpha=0.5)
plt.show()

#fp = '/Users/josh/projects/intro/ngc3621_12m+7m+tp_co21_90pc_props.fits.bz2'
#fp = '/data/kant/0/sun.1608/PHANGS/ALMA/v3p4-CPROPS/STv1p5/90pc_homogenized/ngc3621_12m+7m+tp_co21_90pc_props.fits.bz2'
source = 'ngc3621'
res = np.array([60,90,120,150])
for i in range(len(res)):
    cat = pd.read_csv(str(source)+'_'+str(res[i])+'pc_cloud_stats.csv')
    fp = '/Users/josh/projects/intro/ngc3621_12m+7m+tp_co21_'+str(res[i])+'pc_props.fits.bz2'
    tab = Table.read(fp)
    x = np.array(tab['XCTR_DEG'])
    y = np.array(tab['YCTR_DEG'])
    xs = x*np.cos(np.deg2rad(x))
    vel = tab['VCTR_KMS']

    ### General Plot of Peaks ###
    plt.scatter(xs,y, c=vel, cmap='RdBu', alpha=0.75)
    plt.xlabel('R.A.')
    plt.ylabel('DEC')
    plt.title('NGC3621 ('+str(res[i])+'pc_homogenized) - Dist. Corrected', fontsize=15)
    plt.colorbar(label='Velocity km/s')
    plt.savefig(loc+'/plots/ngc3621_peaks'+str(res[i])+'pc.pdf')
    plt.close()
    print(str(res[i])+'pc peak')

    ### Plotting nearest neighbor separation (PC) ###
    plt.hist(cat.min_dist, bins = (np.linspace(0, 500, 10)), label=r'Eucledian (N = '+str(len(x))+')', alpha = 0.5)
    #14 sources above 500pc
    plt.xlabel('Distance to Nearest Neighbor (pc)')
    plt.ylabel('Counts')
    plt.legend()
    plt.title('NGC3621 ('+str(res[i])+'pc_homogenized)')
    plt.savefig(loc+'/plots/ngc3621_pc_hist'+str(res[i])+'pc.pdf')
    plt.close()
    plt.hist(cat.min_dist2nd, bins = (np.linspace(0, 700, 14)), label=r'Eucledian (N = '+str(len(x))+')', alpha = 0.5)
    plt.xlabel('Distance to 2nd Nearest Neighbor (pc)')
    plt.ylabel('Counts')
    plt.legend()
    plt.title('NGC3621 ('+str(res[i])+'pc_homogenized)')
    plt.savefig(loc+'/plots/ngc3621_pc_hist_2nd'+str(res[i])+'pc.pdf')
    plt.close()
    plt.hist(cat.min_dist3rd, bins = (np.linspace(0, 700, 14)), label=r'Eucledian (N = '+str(len(x))+')', alpha = 0.5)
    plt.xlabel('Distance to 3rd Nearest Neighbor (pc)')
    plt.ylabel('Counts')
    plt.legend()
    plt.title('NGC3621 ('+str(res[i])+'pc_homogenized)')
    plt.savefig(loc+'/plots/ngc3621_pc_hist_3rd'+str(res[i])+'pc.pdf')
    plt.close()
    print(str(res[i])+'pc hist')

    ### Plotting nearest neighbor separation (beam) ###
    plt.hist(cat.beam_sep_nn, bins = (np.linspace(0, 6, 12)), label=r'Eucledian (N = '+str(len(x))+')', alpha = 0.5)
    plt.xlabel('Distance to Nearest Neighbor (beams)')
    plt.ylabel('Counts')
    plt.legend()
    plt.title('NGC3621 ('+str(res[i])+'pc_homogenized)')
    plt.savefig(loc+'/plots/ngc3621_beam_hist'+str(res[i])+'pc.pdf')
    plt.close()
    plt.hist(cat.beam_sep_nn2, bins = (np.linspace(0, 8, 16)), label=r'Eucledian (N = '+str(len(x))+')', alpha = 0.5)
    plt.xlabel('Distance to 2nd Nearest Neighbor (beams)')
    plt.ylabel('Counts')
    plt.legend()
    plt.title('NGC3621 ('+str(res[i])+'pc_homogenized)')
    plt.savefig(loc+'/plots/ngc3621_beam_hist_2nd'+str(res[i])+'pc.pdf')
    plt.close()
    plt.hist(cat.beam_sep_nn3, bins = (np.linspace(0, 10, 20)), label=r'Eucledian (N = '+str(len(x))+')', alpha = 0.5)
    plt.xlabel('Distance to 3rd Nearest Neighbor (beams)')
    plt.ylabel('Counts')
    plt.legend()
    plt.title('NGC3621 ('+str(res[i])+'pc_homogenized)')
    plt.savefig(loc+'/plots/ngc3621_beam_hist_3rd'+str(res[i])+'pc.pdf')
    plt.close()
    print(str(res[i])+'pc beam sep hist')

    ### Plotting nearest neighbor separation (mean cloud distance) ###
    plt.hist(cat.mean_cloud_sep_nn, bins = (np.linspace(0, 15, 15)), label=r'Eucledian (N = '+str(len(x))+')', alpha = 0.5)
    plt.xlabel('Distance to Nearest Neighbor (mean cloud radii)')
    plt.ylabel('Counts')
    plt.legend()
    plt.title('NGC3621 ('+str(res[i])+'pc_homogenized)')
    plt.savefig(loc+'/plots/ngc3621_cloudsize_hist'+str(res[i])+'pc.pdf')
    plt.close()
    plt.hist(cat.mean_cloud_sep_nn2, bins = (np.linspace(0, 15, 15)), label=r'Eucledian (N = '+str(len(x))+')', alpha = 0.5)
    plt.xlabel('Distance to 2nd Nearest Neighbor (mean cloud radii)')
    plt.ylabel('Counts')
    plt.legend()
    plt.title('NGC3621 ('+str(res[i])+'pc_homogenized)')
    plt.savefig(loc+'/plots/ngc3621_cloudsize_hist_2nd'+str(res[i])+'pc.pdf')
    plt.close()
    plt.hist(cat.mean_cloud_sep_nn3, bins = (np.linspace(0, 15, 15)), label=r'Eucledian (N = '+str(len(x))+')', alpha = 0.5)
    plt.xlabel('Distance to 3rd Nearest Neighbor (mean cloud radii)')
    plt.ylabel('Counts')
    plt.legend()
    plt.title('NGC3621 ('+str(res[i])+'pc_homogenized)')
    plt.savefig(loc+'/plots/ngc3621_cloudsize_hist_3rd'+str(res[i])+'pc.pdf')
    plt.close()
    print(str(res[i])+'pc cloud size sep hist')
