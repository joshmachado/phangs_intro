from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import pandas as pd
import matplotlib.lines as mlines

loc = '/Users/josh/projects/intro/'
source = 'ngc2835'

sep = np.zeros([3,3,4,2])
msep = np.zeros([3,3,4,2])
err = np.zeros([3,3,2,4])
merr = np.zeros([3,3,2,4])
res = [60,90,120,150]
for i in range(len(res)):
    cat = pd.read_csv(loc+source+'/'+source+'_'+str(res[i])+'pc_cloud_stats.csv')
    match_cat = pd.read_csv(loc+source+'/matched/'+str(source)+'_'+str(res[i])+'pc_cloud_stats.csv')

    #Separation in PC
    pc = cat['min_dist']
    pc_matched = match_cat['min_dist']
    pc2 = cat['min_dist2nd']
    pc2_matched = match_cat['min_dist2nd']
    pc3 = cat['min_dist3rd']
    pc3_matched = match_cat['min_dist3rd']
    sep[0,0,i] = [res[i], np.percentile(pc, 50)]
    msep[0,0,i] = [res[i], np.percentile(pc_matched, 50)]
    err[0,0,:,i] = [np.percentile(pc, 16), np.percentile(pc, 84)]
    merr[0,0,:,i] = [np.percentile(pc_matched, 16), np.percentile(pc_matched, 84)]
    sep[0,1,i] = [res[i], np.percentile(pc2, 50)]
    msep[0,1,i] = [res[i], np.percentile(pc2_matched, 50)]
    err[0,1,:,i] = [np.percentile(pc2, 16), np.percentile(pc2, 84)]
    merr[0,1,:,i] = [np.percentile(pc2_matched, 16), np.percentile(pc2_matched, 84)]
    sep[0,2,i] = [res[i], np.percentile(pc3, 50)]
    msep[0,2,i] = [res[i], np.percentile(pc3_matched, 50)]
    err[0,2,:,i] = [np.percentile(pc3, 16), np.percentile(pc3, 84)]
    merr[0,2,:,i] = [np.percentile(pc3_matched, 16), np.percentile(pc3_matched, 84)]

    #Separation in Beams
    beam = cat['beam_sep_nn']
    beam_matched = match_cat['beam_sep_nn']
    beam2 = cat['beam_sep_nn2']
    beam2_matched = match_cat['beam_sep_nn2']
    beam3 = cat['beam_sep_nn3']
    beam3_matched = match_cat['beam_sep_nn3']
    sep[1,0,i] = [res[i], np.percentile(beam, 50)]
    msep[1,0,i] = [res[i], np.percentile(beam_matched, 50)]
    err[1,0,:,i] = [np.percentile(beam, 16), np.percentile(beam, 84)]
    merr[1,0,:,i] = [np.percentile(beam_matched, 16), np.percentile(beam_matched, 84)]
    sep[1,1,i] = [res[i], np.percentile(beam2, 50)]
    msep[1,1,i] = [res[i], np.percentile(beam2_matched, 50)]
    err[1,1,:,i] = [np.percentile(beam2, 16), np.percentile(beam2, 84)]
    merr[1,1,:,i] = [np.percentile(beam2_matched, 16), np.percentile(beam2_matched, 84)]
    sep[1,2,i] = [res[i], np.percentile(beam3, 50)]
    msep[1,2,i] = [res[i], np.percentile(beam3_matched, 50)]
    err[1,2,:,i] = [np.percentile(beam3, 16), np.percentile(beam3, 84)]
    merr[1,2,:,i] = [np.percentile(beam3_matched, 16), np.percentile(beam3_matched, 84)]

    #Separation in Mean Cloud Radii
    cloud = cat['mean_cloud_sep_nn']
    cloud_matched = match_cat['mean_cloud_sep_nn']
    cloud2 = cat['mean_cloud_sep_nn2']
    cloud2_matched = match_cat['mean_cloud_sep_nn2']
    cloud3 = cat['mean_cloud_sep_nn3']
    cloud3_matched = match_cat['mean_cloud_sep_nn3']
    sep[2,0,i] = [res[i], np.percentile(cloud, 50)]
    msep[2,0,i] = [res[i], np.percentile(cloud_matched, 50)]
    err[2,0,:,i] = [np.percentile(cloud, 16), np.percentile(cloud, 84)]
    merr[2,0,:,i] = [np.percentile(cloud_matched, 16), np.percentile(cloud_matched, 84)]
    sep[2,1,i] = [res[i], np.percentile(cloud2, 50)]
    msep[2,1,i] = [res[i], np.percentile(cloud2_matched, 50)]
    err[2,1,:,i] = [np.percentile(cloud2, 16), np.percentile(cloud2, 84)]
    merr[2,1,:,i] = [np.percentile(cloud2_matched, 16), np.percentile(cloud2_matched, 84)]
    sep[2,2,i] = [res[i], np.percentile(cloud3, 50)]
    msep[2,2,i] = [res[i], np.percentile(cloud3_matched, 50)]
    err[2,2,:,i] = [np.percentile(cloud3, 16), np.percentile(cloud3, 84)]
    merr[2,2,:,i] = [np.percentile(cloud3_matched, 16), np.percentile(cloud3_matched, 84)]


#FORMATTING LEGEND
first = mlines.Line2D([], [], color='#377eb8', lw=1, label='1st Nearest Neighbor')
second = mlines.Line2D([], [], color='#4daf4a', lw=1, label='2nd Nearest Neighbor')
third = mlines.Line2D([], [], color='#a65628', lw=1, label='3rd Nearest Neighbor')
homogenized = mlines.Line2D([], [], color='black', marker='o', lw=0,  label='Homogenized')
matched = mlines.Line2D([], [], color='black', marker='s', lw=0,  label='Matched')

#PLOTTING PC SEPERATION DATA
plt.figure(figsize=(10,10))
plt.scatter(sep[0,0,:,0]-8, sep[0,0,:,1], c='#377eb8', label='1st Homogenized')
plt.errorbar(sep[0,0,:,0]-8, sep[0,0,:,1], c='#377eb8', xerr=None, yerr=err[0,0], capsize=4, ls='none')
plt.scatter(msep[0,0,:,0]-6, msep[0,0,:,1], c='#377eb8', marker='s', label='1st Matched')
plt.errorbar(msep[0,0,:,0]-6, msep[0,0,:,1], c='#377eb8', xerr=None, yerr=merr[0,0], capsize=4, ls='none')
plt.scatter(sep[0,1,:,0]-1, sep[0,1,:,1], c='#4daf4a', label='2nd Homogenized')
plt.errorbar(sep[0,1,:,0]-1, sep[0,1,:,1], c='#4daf4a', xerr=None, yerr=err[0,1], capsize=4, ls='none')
plt.scatter(msep[0,1,:,0]+1, msep[0,1,:,1], c='#4daf4a', marker='s', label='2nd Matched')
plt.errorbar(msep[0,1,:,0]+1, msep[0,1,:,1], c='#4daf4a', xerr=None, yerr=merr[0,1], capsize=4, ls='none')
plt.scatter(sep[0,2,:,0]+5, sep[0,2,:,1], c='#a65628', label='3rd Homogenized')
plt.errorbar(sep[0,2,:,0]+5, sep[0,2,:,1], c='#a65628', xerr=None, yerr=err[0,2], capsize=4, ls='none')
plt.scatter(msep[0,2,:,0]+7, msep[0,2,:,1], c='#a65628', marker='s', label='3rd Matched')
plt.errorbar(msep[0,2,:,0]+7, msep[0,2,:,1], c='#a65628', xerr=None, yerr=merr[0,2], capsize=4, ls='none')

#AXIS & LABELS
plt.legend(handles=[first, second, third, homogenized, matched])
plt.xlabel('Resolution (pc)')
plt.ylabel('Distance to Nearest Neighbor (pc)')
plt.xticks([60,90,120,150])
plt.title(str(source).upper(), fontsize=16)
plt.savefig(loc+source+'/plots/'+str(source)+'_pc_sep_summary.pdf')
plt.close()

#PLOTTING BEAM SEPERATION DATA
plt.figure(figsize=(10,10))
plt.scatter(sep[1,0,:,0]-8, sep[1,0,:,1], c='#377eb8', label='1st Homogenized')
plt.errorbar(sep[1,0,:,0]-8, sep[1,0,:,1], c='#377eb8', xerr=None, yerr=err[1,0], capsize=4, ls='none')
plt.scatter(msep[1,0,:,0]-6, msep[1,0,:,1], c='#377eb8', marker='s', label='1st Matched')
plt.errorbar(msep[1,0,:,0]-6, msep[1,0,:,1], c='#377eb8', xerr=None, yerr=merr[1,0], capsize=4, ls='none')
plt.scatter(sep[1,1,:,0]-1, sep[1,1,:,1], c='#4daf4a', label='2nd Homogenized')
plt.errorbar(sep[1,1,:,0]-1, sep[1,1,:,1], c='#4daf4a', xerr=None, yerr=err[1,1], capsize=4, ls='none')
plt.scatter(msep[1,1,:,0]+1, msep[1,1,:,1], c='#4daf4a', marker='s', label='2nd Matched')
plt.errorbar(msep[1,1,:,0]+1, msep[1,1,:,1], c='#4daf4a', xerr=None, yerr=merr[1,1], capsize=4, ls='none')
plt.scatter(sep[1,2,:,0]+5, sep[1,2,:,1], c='#a65628', label='3rd Homogenized')
plt.errorbar(sep[1,2,:,0]+5, sep[1,2,:,1], c='#a65628', xerr=None, yerr=err[1,2], capsize=4, ls='none')
plt.scatter(msep[1,2,:,0]+7, msep[1,2,:,1], c='#a65628', marker='s', label='3rd Matched')
plt.errorbar(msep[1,2,:,0]+7, msep[1,2,:,1], c='#a65628', xerr=None, yerr=merr[1,2], capsize=4, ls='none')

#AXIS & LABELS
plt.legend(handles=[first, second, third, homogenized, matched])
plt.xlabel('Resolution (pc)')
plt.ylabel('Distance to Nearest Neighbor (beams)')
#plt.ylim(0,np.max(merr[1,2]+5))
plt.xticks([60,90,120,150])
plt.title(str(source).upper(), fontsize=16)
plt.savefig(loc+source+'/plots/'+str(source)+'_beam_sep_summary.pdf')
plt.close()

#PLOTTING MEAN CLOUD RADII SEPERATION DATA
plt.figure(figsize=(10,10))
plt.scatter(sep[2,0,:,0]-8, sep[2,0,:,1], c='#377eb8', label='1st Homogenized')
plt.errorbar(sep[2,0,:,0]-8, sep[2,0,:,1], c='#377eb8', xerr=None, yerr=err[2,0], capsize=4, ls='none')
plt.scatter(msep[2,0,:,0]-6, msep[2,0,:,1], c='#377eb8', marker='s', label='1st Matched')
plt.errorbar(msep[2,0,:,0]-6, msep[2,0,:,1], c='#377eb8', xerr=None, yerr=merr[2,0], capsize=4, ls='none')
plt.scatter(sep[2,1,:,0]-1, sep[2,1,:,1], c='#4daf4a', label='2nd Homogenized')
plt.errorbar(sep[2,1,:,0]-1, sep[2,1,:,1], c='#4daf4a', xerr=None, yerr=err[2,1], capsize=4, ls='none')
plt.scatter(msep[2,1,:,0]+1, msep[2,1,:,1], c='#4daf4a', marker='s', label='2nd Matched')
plt.errorbar(msep[2,1,:,0]+1, msep[2,1,:,1], c='#4daf4a', xerr=None, yerr=merr[2,1], capsize=4, ls='none')
plt.scatter(sep[2,2,:,0]+5, sep[2,2,:,1], c='#a65628', label='3rd Homogenized')
plt.errorbar(sep[2,2,:,0]+5, sep[2,2,:,1], c='#a65628', xerr=None, yerr=err[2,2], capsize=4, ls='none')
plt.scatter(msep[2,2,:,0]+7, msep[2,2,:,1], c='#a65628', marker='s', label='3rd Matched')
plt.errorbar(msep[2,2,:,0]+7, msep[2,2,:,1], c='#a65628', xerr=None, yerr=merr[2,2], capsize=4, ls='none')

#AXIS & LABELS
plt.legend(handles=[first, second, third, homogenized, matched])
plt.xlabel('Resolution (pc)')
plt.ylabel('Distance to Nearest Neighbor (mean cloud radii)')
plt.xticks([60,90,120,150])
plt.title(str(source).upper(), fontsize=16)
plt.savefig(loc+source+'/plots/'+str(source)+'_cloud_sep_summary.pdf')
plt.close()
