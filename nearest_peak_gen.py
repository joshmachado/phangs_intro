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
    fp = '/Users/josh/projects/intro/ngc3621_12m+7m+tp_co21_'+str(res[i])+'pc_props.fits.bz2'
    tab = Table.read(fp)
    cloudnum = np.array(tab['CLOUDNUM'])
    distance = np.array(tab['DISTANCE_PC'])
    x = np.array(tab['XCTR_DEG'])
    y = np.array(tab['YCTR_DEG'])
    xs = x*np.cos(np.deg2rad(x)) #Spherical correction? RA = RAcos(Dec)
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
        ind3 = np.partition(dist[k], 2)[2] #index of 3rd nearest neighbor
        mindist[k] = dist[k,int(ind[0])]
        nn[k] = int(ind[0])+1
        mindist2[k] = dist[k,int(ind2[0])]
        nn2[k] = int(ind2[0])+1
        mindist3[k] = dist[k,int(ind3[0])]
        nn3[k] = int(ind3[0])+1
    for k in range(len(x)):
        ind = np.where(corr_dist[k] == np.nanmin(corr_dist[k]))
        mindist_corr[k] = corr_dist[k,int(ind[0])]
        nn_corr[k] = int(ind[0])+1
        #ind2 = np.where(corr_dist[k])

    true_dist = np.deg2rad(mindist)*distance[0]
    true_dist_corr = np.deg2rad(mindist_corr)*distance[0]
    cat = pd.DataFrame({'cloudnum':cloudnum, 'x':x, 'corr_x':xs, 'y':y, 'nearest_neighbor':nn,
    'nearest_neighbor_corr':nn_corr, 'min_distance':mindist, 'min_distance_corr':mindist_corr, 'true_dist':true_dist, 'true_dist_corr':true_dist_corr})
    cat.to_csv('ngc3621_'+str(res[i])+'pc.csv')
    print(str(res[i])+'pc done')
