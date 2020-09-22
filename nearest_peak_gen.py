from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import pandas as pd

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
    ys = y*np.cos(x) #Spherical correction? RA = RAcos(Dec)
    peaks = np.array([cloudnum,x,y])
    dist = np.zeros((len(x),len(x)))
    corr_dist = np.zeros((len(x),len(x)))
    min = np.zeros(len(x)) #stores index of nearest neighbor
    min_corr = np.zeros(len(x))
    nn = np.zeros(len(x)) #stores CloudNum of nearest neighbor
    nn_corr = np.zeros(len(x))
    mindist = np.zeros(len(x)) #stores distance to nearest neighbor
    mindist_corr = np.zeros(len(x))
    ## row: peaks[:,0] returns cloudnum, x,y
    k = 0
    j = 0
    while k < len(x):
        for j in range(len(x)):
            dist[k,j] = np.sqrt(np.square(x[k] - x[j]) + np.square(y[k] - y[j]))
            corr_dist[k,j] = np.sqrt(np.square(x[k] - x[j]) + np.square(ys[k] - ys[j]))
        k+=1
    dist[dist == 0] = np.nan
    corr_dist[corr_dist == 0] = np.nan
    for k in range(len(x)):
        ind = np.where(dist[k] == np.nanmin(dist[k]))
        mindist[k] = dist[k,int(ind[0])]
        nn[k] = int(ind[0])+1
    for k in range(len(x)):
        ind = np.where(corr_dist[k] == np.nanmin(corr_dist[k]))
        mindist_corr[k] = corr_dist[k,int(ind[0])]
        nn_corr[k] = int(ind[0])+1

    true_dist = mindist*distance[0]
    true_dist_corr = mindist_corr*distance[0]
    cat = pd.DataFrame({'cloudnum':cloudnum, 'x':x, 'y':y, 'corrected_y':ys, 'nearest_neighbor':nn,
    'nearest_neighbor_corr':nn_corr, 'min_distance':mindist, 'min_distance_corr':mindist_corr, 'true_dist':true_dist, 'true_dist_corr':true_dist_corr})
    cat.to_csv('ngc3621_'+str(res[i])+'pc.csv')
