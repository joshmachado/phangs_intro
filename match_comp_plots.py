from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import pandas as pd
loc = '/Users/josh/projects/intro'
match_loc = '/Users/josh/projects/intro/matched'

#fp = '/Users/josh/projects/intro/ngc3621_12m+7m+tp_co21_90pc_props.fits.bz2'
#fp = '/data/kant/0/sun.1608/PHANGS/ALMA/v3p4-CPROPS/STv1p5/90pc_homogenized/ngc3621_12m+7m+tp_co21_90pc_props.fits.bz2'
res = np.array([60,90,120,150])
for i in range(len(res)):
    cat = pd.read_csv('ngc3621_'+str(res[i])+'pc_cloud_stats.csv')
    fp = '/Users/josh/projects/intro/ngc3621_12m+7m+tp_co21_'+str(res[i])+'pc_props.fits.bz2'
    tab = Table.read(fp)
    x = np.array(tab['XCTR_DEG'])
    y = np.array(tab['YCTR_DEG'])
    xs = x*np.cos(np.deg2rad(x))
    vel = tab['VCTR_KMS']

    ##MATCHED
    match_cat = pd.read_csv(match_loc+'/ngc3621_'+str(res[i])+'pc_cloud_stats.csv')
    match_fp = '/Users/josh/projects/intro/matched/ngc3621_12m+7m+tp_co21_'+str(res[i])+'pc_props.fits.bz2'
    match_tab = Table.read(match_fp)
    match_x = np.array(match_tab['XCTR_DEG'])
    match_y = np.array(match_tab['YCTR_DEG'])
    match_xs = match_x*np.cos(np.deg2rad(match_x))
    match_vel = match_tab['VCTR_KMS']


    plt.scatter(xs,y, alpha=0.75, label='Homogenized')
    plt.scatter(match_xs, match_y, alpha = 0.25, label='Matched')
    plt.xlabel('R.A.')
    plt.ylabel('DEC')
    plt.title('NGC3621 ('+str(res[i])+'pc Matched vs. Homogenized)', fontsize=15)
    plt.legend()
    plt.savefig(match_loc+'/plots/ngc3621_comp_peaks'+str(res[i])+'pc.pdf')
    plt.close()
    print(str(res[i])+'pc peak')

    plt.hist(cat.min_dist, bins = (np.linspace(0, 700, 14)), label=r'Homogenized (N = '+str(len(x))+')', alpha = 0.5)
    #14 sources above 500pc
    plt.hist(match_cat.min_dist, bins = (np.linspace(0, 700, 14)), label=r'Matched (N = '+str(len(match_x))+')', alpha = 0.5)
    plt.xlabel('Distance to Nearest Neighbor (pc)')
    plt.ylabel('Counts')
    plt.legend()
    plt.title('NGC3621 ('+str(res[i])+'pc)')
    plt.savefig(match_loc+'/plots/ngc3621_hist'+str(res[i])+'pc.pdf')
    plt.close()

    plt.hist(cat.min_dist2nd, bins = (np.linspace(0, 700, 14)), label=r'Homogenized (N = '+str(len(x))+')', alpha = 0.5)
    plt.hist(match_cat.min_dist2nd, bins = (np.linspace(0, 700, 14)), label=r'Matched (N = '+str(len(match_x))+')', alpha = 0.5)
    plt.xlabel('Distance to 2nd Nearest Neighbor (pc)')
    plt.ylabel('Counts')
    plt.legend()
    plt.title('NGC3621 ('+str(res[i])+'pc)')
    plt.savefig(match_loc+'/plots/ngc3621_hist_2nd'+str(res[i])+'pc.pdf')
    plt.close()

    plt.hist(cat.min_dist3rd, bins = (np.linspace(0, 700, 14)), label=r'Homogenized (N = '+str(len(x))+')', alpha = 0.5)
    plt.hist(match_cat.min_dist3rd, bins = (np.linspace(0, 700, 14)), label=r'Matched (N = '+str(len(match_x))+')', alpha = 0.5)
    plt.xlabel('Distance to 3rd Nearest Neighbor (pc)')
    plt.ylabel('Counts')
    plt.legend()
    plt.title('NGC3621 ('+str(res[i])+'pc)')
    plt.savefig(match_loc+'/plots/ngc3621_hist_3rd'+str(res[i])+'pc.pdf')
    plt.close()

    ##plot cloud sep and beam separation
    plt.hist(cat.beam_sep_nn, bins = (np.linspace(0, 5, 10)), label=r'Homogenized (N = '+str(len(x))+')', alpha = 0.5)
    plt.hist(match_cat.beam_sep_nn, bins = (np.linspace(0, 5, 10)), label=r'Matched (N = '+str(len(match_x))+')', alpha = 0.5)
    plt.xlabel('Beam Separation to Nearest Neighbor')
    plt.ylabel('Counts')
    plt.legend()
    plt.title('NGC3621 ('+str(res[i])+'pc)')
    plt.savefig(match_loc+'/plots/ngc3621_beamsep'+str(res[i])+'pc.pdf')
    plt.close()

    plt.hist(cat.beam_sep_nn2, bins = (np.linspace(0, 7, 14)), label=r'Homogenized (N = '+str(len(x))+')', alpha = 0.5)
    plt.hist(match_cat.beam_sep_nn2, bins = (np.linspace(0, 7, 14)), label=r'Matched (N = '+str(len(match_x))+')', alpha = 0.5)
    plt.xlabel('Beam Separation to 2nd Nearest Neighbor')
    plt.ylabel('Counts')
    plt.legend()
    plt.title('NGC3621 ('+str(res[i])+'pc)')
    plt.savefig(match_loc+'/plots/ngc3621_beamsep2nd'+str(res[i])+'pc.pdf')
    plt.close()

    plt.hist(cat.beam_sep_nn3, bins = (np.linspace(0, 10, 20)), label=r'Homogenized (N = '+str(len(x))+')', alpha = 0.5)
    plt.hist(match_cat.beam_sep_nn3, bins = (np.linspace(0, 10, 20)), label=r'Matched (N = '+str(len(match_x))+')', alpha = 0.5)
    plt.xlabel('Beam Separation to 3rd Nearest Neighbor')
    plt.ylabel('Counts')
    plt.legend()
    plt.title('NGC3621 ('+str(res[i])+'pc)')
    plt.savefig(match_loc+'/plots/ngc3621_beamsep3rd'+str(res[i])+'pc.pdf')
    plt.close()

    plt.hist(cat.mean_cloud_sep_nn, bins = (np.linspace(0, 7, 14)), label=r'Homogenized (N = '+str(len(x))+')', alpha = 0.5)
    plt.hist(match_cat.mean_cloud_sep_nn, bins = (np.linspace(0, 7, 14)), label=r'Matched (N = '+str(len(match_x))+')', alpha = 0.5)
    plt.xlabel('Mean Cloud Radii Separation to Nearest Neighbor')
    plt.ylabel('Counts')
    plt.legend()
    plt.title('NGC3621 ('+str(res[i])+'pc)')
    plt.savefig(match_loc+'/plots/ngc3621_cloudsep'+str(res[i])+'pc.pdf')
    plt.close()

    plt.hist(cat.mean_cloud_sep_nn2, bins = (np.linspace(0, 7, 14)), label=r'Homogenized (N = '+str(len(x))+')', alpha = 0.5)
    plt.hist(match_cat.mean_cloud_sep_nn2, bins = (np.linspace(0, 7, 14)), label=r'Matched (N = '+str(len(match_x))+')', alpha = 0.5)
    plt.xlabel('Mean Cloud Radii Separation to 2nd Nearest Neighbor')
    plt.ylabel('Counts')
    plt.legend()
    plt.title('NGC3621 ('+str(res[i])+'pc)')
    plt.savefig(match_loc+'/plots/ngc3621_cloudsep2nd'+str(res[i])+'pc.pdf')
    plt.close()

    plt.hist(cat.mean_cloud_sep_nn3, bins = (np.linspace(0, 10, 20)), label=r'Homogenized (N = '+str(len(x))+')', alpha = 0.5)
    plt.hist(match_cat.mean_cloud_sep_nn3, bins = (np.linspace(0, 10, 20)), label=r'Matched (N = '+str(len(match_x))+')', alpha = 0.5)
    plt.xlabel('Mean Cloud Radii Separation to 3rd Nearest Neighbor')
    plt.ylabel('Counts')
    plt.legend()
    plt.title('NGC3621 ('+str(res[i])+'pc)')
    plt.savefig(match_loc+'/plots/ngc3621_cloudsep3rd'+str(res[i])+'pc.pdf')
    plt.close()
    print(str(res[i])+'pc hist')
