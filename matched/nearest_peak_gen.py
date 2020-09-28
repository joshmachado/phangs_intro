from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import pandas as pd

#cp /data/tycho/0/sun.1608/PHANGS/ALMA/v3p4-CPROPS/STv1p5/150pc_matched/ngc3621_12m+7m+tp_co21_150pc_props.fits.bz2 .

loc = '/Users/josh/projects/intro/matched/'
#fp = '/Users/josh/projects/intro/ngc3621_12m+7m+tp_co21_90pc_props.fits.bz2'
#fp = '/data/kant/0/sun.1608/PHANGS/ALMA/v3p4-CPROPS/STv1p5/90pc_homogenized/ngc3621_12m+7m+tp_co21_90pc_props.fits.bz2'
res = np.array([60,90,120,150])
mean_dist = np.zeros([4,3])
beam_sep = np.zeros([4,3])
percentiles = np.zeros([4,6]) #90th & 75th
for i in range(len(res)):
    fp = '/Users/josh/projects/intro/matched/ngc3621_12m+7m+tp_co21_'+str(res[i])+'pc_props.fits.bz2'
    tab = Table.read(fp)
    cloudnum = np.array(tab['CLOUDNUM'])
    distance = np.array(tab['DISTANCE_PC'])
    x = np.array(tab['XCTR_DEG'])
    y = np.array(tab['YCTR_DEG'])
    xs = x*np.cos(np.deg2rad(y)) #Spherical correction? RA = RAcos(Dec)
    peaks = np.array([cloudnum,x,y])
    dist = np.zeros((len(x),len(x)))
    corr_dist = np.zeros((len(x),len(x)))
    min = np.zeros(len(x)) #stores index of nearest neighbor
    min_corr = np.zeros(len(x))
    nn = np.zeros(len(x)) #stores CloudNum of nearest neighbor
    nn_corr = np.zeros(len(x))
    mindist = np.zeros(len(x)) #stores distance to nearest neighbor
    mindist_corr = np.zeros(len(x))
    nn2 = np.zeros(len(x)) #stores CloudNum of 2nd nearest neighbor
    nn2_corr = np.zeros(len(x))
    mindist2 = np.zeros(len(x)) #stores distance to 2nd nearest neighbor
    mindist2_corr = np.zeros(len(x))
    nn3 = np.zeros(len(x)) #stores CloudNum of 3rd nearest neighbor
    nn3_corr = np.zeros(len(x))
    mindist3 = np.zeros(len(x)) #stores distance to 3rd nearest neighbor
    mindist3_corr = np.zeros(len(x))
    ## row: peaks[:,0] returns cloudnum, x,y
    k = 0
    j = 0
    while k < len(x):
        for j in range(len(x)):
            dist[k,j] = np.sqrt(np.square(x[k] - x[j]) + np.square(y[k] - y[j]))
            corr_dist[k,j] = np.sqrt(np.square(xs[k] - xs[j]) + np.square(y[k] - y[j]))
        k+=1
    dist[dist == 0] = np.nan
    corr_dist[corr_dist == 0] = np.nan
    for k in range(len(x)):
        ind = np.where(dist[k] == np.nanmin(dist[k]))
        ind2 = np.where(dist[k] == np.partition(dist[k], 1)[1]) #index of 2nd nearest neighbor
        ind3 = np.where(dist[k] == np.partition(dist[k], 2)[2]) #index of 3rd nearest neighbor
        mindist[k] = np.deg2rad(dist[k,int(ind[0])])*distance[0]
        nn[k] = int(ind[0])+1
        mindist2[k] = np.deg2rad(dist[k,int(ind2[0])])*distance[0]
        nn2[k] = int(ind2[0])+1
        mindist3[k] = np.deg2rad(dist[k,int(ind3[0])])*distance[0]
        nn3[k] = int(ind3[0])+1
    for k in range(len(x)):
        ind = np.where(corr_dist[k] == np.nanmin(corr_dist[k]))
        mindist_corr[k] = np.deg2rad(corr_dist[k,int(ind[0])])*distance[0]
        nn_corr[k] = int(ind[0])+1
        ind2 = np.where(corr_dist[k] == np.partition(corr_dist[k], 1)[1])
        mindist2_corr[k] = np.deg2rad(corr_dist[k,int(ind2[0])])*distance[0]
        nn2_corr[k] = int(ind2[0])+1
        ind3 = np.where(corr_dist[k] == np.partition(corr_dist[k], 2)[2])
        mindist3_corr[k] = np.deg2rad(corr_dist[k,int(ind3[0])])*distance[0]
        nn3_corr[k] = int(ind3[0])+1

    true_dist = np.deg2rad(mindist)*distance[0]
    true_dist_corr = np.deg2rad(mindist_corr)*distance[0]
    true_dist2 = np.deg2rad(mindist2)*distance[0]
    true_dist2_corr = np.deg2rad(mindist2_corr)*distance[0]
    true_dist3 = np.deg2rad(mindist3)*distance[0]
    true_dist3_corr = np.deg2rad(mindist3_corr)*distance[0]
    cat = pd.DataFrame({'cloudnum':cloudnum, 'x':x, 'corr_x':xs, 'y':y, 'nearest_neighbor':nn,
    'nearest_neighbor_corr':nn_corr, 'min_distance':mindist, 'min_distance_corr':mindist_corr, 'true_dist':true_dist, 'true_dist_corr':true_dist_corr,
    'nearest_neighbor2':nn2, 'nearest_neighbor2_corr':nn2_corr, 'min_distance2':mindist2, 'min_distance_corr2':mindist2_corr, 'true_dist2_corr':true_dist2_corr,
    'nearest_neighbor3':nn3, 'nearest_neighbor3_corr':nn3_corr, 'min_distance3':mindist3, 'min_distance_corr3':mindist3_corr, 'true_dist3_corr':true_dist3_corr})
    cat.to_csv(loc+'ngc3621_'+str(res[i])+'pc.csv', index=False)

    ###LETS GATHER SOME STATS###
    mean_dist[i] = [np.mean(mindist_corr), np.mean(mindist2_corr), np.mean(mindist3_corr)]
    percentiles[i] = [np.percentile(true_dist_corr, 90),np.percentile(true_dist2_corr, 90),np.percentile(true_dist3_corr, 90),
    np.percentile(true_dist_corr, 75),np.percentile(true_dist2_corr, 75),np.percentile(true_dist3_corr, 75)]
    beam_sep[i] = [np.mean(mindist_corr)/tab['BEAMFWHM_PC'][0], np.mean(mindist2_corr)/tab['BEAMFWHM_PC'][0], np.mean(mindist3_corr)/tab['BEAMFWHM_PC'][0]]
    print(str(res[i])+'pc done')
print('Gathering stats...')

peak_stats = pd.DataFrame({'res_pc':res, 'mean_dist_nn':mean_dist[:,0], 'mean_dist_nn2':mean_dist[:,1],
'mean_dist_nn3':mean_dist[:,2], '90th_nn':percentiles[:,0], '90th_nn2':percentiles[:,1], '90th_nn3':percentiles[:,2],
'75th_nn':percentiles[:,3], '75th_nn2':percentiles[:,4], '75th_nn3':percentiles[:,5], 'mean_beam_sep':beam_sep[:,0],
'mean_beam_sep2':beam_sep[:,1], 'mean_beam_sep3':beam_sep[:,2]})

peak_stats.to_csv(loc+'peak_stats.csv', index=False)
