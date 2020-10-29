### Comparing Max vs. Centers ###
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import pandas as pd
import matplotlib.lines as mlines
import matplotlib.cm as cm


#### NGC3621 ####
#Center - Homogenized
ngc3621 = pd.read_csv('/Users/josh/projects/intro/ngc3621/ngc3621_90pc_cloud_stats.csv')
ngc3621_beams_nn = ngc3621['beam_sep_nn']
ngc3621_beams_nn2 = ngc3621['beam_sep_nn2']
ngc3621_beams_nn3 = ngc3621['beam_sep_nn3']

#Center - Matched
ngc3621_match = pd.read_csv('/Users/josh/projects/intro/ngc3621/matched/ngc3621_90pc_cloud_stats.csv')
ngc3621_match_beams_nn = ngc3621_match['beam_sep_nn']
ngc3621_match_beams_nn2 = ngc3621_match['beam_sep_nn2']
ngc3621_match_beams_nn3 = ngc3621_match['beam_sep_nn3']

#Max
ngc3621_m = pd.read_csv('/Users/josh/projects/intro/ngc3621/max/ngc3621_90pc_cloud_stats.csv')
ngc3621_m_beams_nn = ngc3621_m['beam_sep_nn']
ngc3621_m_beams_nn2 = ngc3621_m['beam_sep_nn2']
ngc3621_m_beams_nn3 = ngc3621_m['beam_sep_nn3']

#### NGC6300 ####
#Center - Homogenized
ngc6300 = pd.read_csv('/Users/josh/projects/intro/ngc6300/ngc6300_90pc_cloud_stats.csv')
ngc6300_beams_nn = ngc6300['beam_sep_nn']
ngc6300_beams_nn2 = ngc6300['beam_sep_nn2']
ngc6300_beams_nn3 = ngc6300['beam_sep_nn3']

#Center - Matched
ngc6300_match = pd.read_csv('/Users/josh/projects/intro/ngc6300/matched/ngc6300_90pc_cloud_stats.csv')
ngc6300_match_beams_nn = ngc6300_match['beam_sep_nn']
ngc6300_match_beams_nn2 = ngc6300_match['beam_sep_nn2']
ngc6300_match_beams_nn3 = ngc6300_match['beam_sep_nn3']

#Max
ngc6300_m = pd.read_csv('/Users/josh/projects/intro/ngc6300/max/ngc6300_90pc_cloud_stats.csv')
ngc6300_m_beams_nn = ngc6300_m['beam_sep_nn']
ngc6300_m_beams_nn2 = ngc6300_m['beam_sep_nn2']
ngc6300_m_beams_nn3 = ngc6300_m['beam_sep_nn3']

#FORMATTING LEGEND
first = mlines.Line2D([], [], color='#377eb8', lw=1, label='1st Nearest Neighbor')
second = mlines.Line2D([], [], color='#4daf4a', lw=1, label='2nd Nearest Neighbor')
third = mlines.Line2D([], [], color='#a65628', lw=1, label='3rd Nearest Neighbor')
centers = mlines.Line2D([], [], color='black', marker='o', lw=0,  label='Center - Homogenized')
matched = mlines.Line2D([], [], color='black', marker='x', lw=0, label='Center - Matched')
maxes = mlines.Line2D([], [], color='black', marker='s', lw=0,  label='Max')

#PLOTTING BEAM SEPERATION DATA

#Centers
plt.figure(figsize=(10,10))
plt.scatter(90.5,  np.median(ngc3621_beams_nn), c='#377eb8', label='1st Nearest Neighbor')
plt.errorbar(90.5, [np.percentile(ngc3621_beams_nn, 16),np.percentile(ngc3621_beams_nn, 84)], c='#377eb8', capsize=4)
plt.scatter(91,  np.median(ngc3621_beams_nn2), c='#4daf4a', label='2nd Nearest Neighbor')
plt.errorbar(91, [np.percentile(ngc3621_beams_nn2, 16),np.percentile(ngc3621_beams_nn2, 84)], c='#4daf4a', capsize=4)
plt.scatter(91.5,  np.median(ngc3621_beams_nn3), c='#a65628', label='3rd Nearest Neighbor')
plt.errorbar(91.5, [np.percentile(ngc3621_beams_nn3, 16),np.percentile(ngc3621_beams_nn3, 84)], c='#a65628', capsize=4)

plt.scatter(90.75,  np.median(ngc3621_match_beams_nn), c='#377eb8', marker='x', label='1st Nearest Neighbor')
plt.errorbar(90.75, [np.percentile(ngc3621_match_beams_nn, 16),np.percentile(ngc3621_match_beams_nn, 84)], c='#377eb8', capsize=4)
plt.scatter(91.2,  np.median(ngc3621_match_beams_nn2), c='#4daf4a', marker='x', label='2nd Nearest Neighbor')
plt.errorbar(91.2, [np.percentile(ngc3621_match_beams_nn2, 16),np.percentile(ngc3621_match_beams_nn2, 84)], c='#4daf4a', capsize=4)
plt.scatter(91.7,  np.median(ngc3621_match_beams_nn3), c='#a65628', marker='x', label='3rd Nearest Neighbor')
plt.errorbar(91.7, [np.percentile(ngc3621_match_beams_nn3, 16),np.percentile(ngc3621_match_beams_nn3, 84)], c='#a65628', capsize=4)

plt.scatter(87.5,  np.median(ngc6300_beams_nn), c='#377eb8', label='1st Nearest Neighbor')
plt.errorbar(87.5, [np.percentile(ngc6300_beams_nn, 16),np.percentile(ngc6300_beams_nn, 84)], c='#377eb8', capsize=4)
plt.scatter(88,  np.median(ngc6300_beams_nn2), c='#4daf4a', label='2nd Nearest Neighbor')
plt.errorbar(88, [np.percentile(ngc6300_beams_nn2, 16),np.percentile(ngc6300_beams_nn2, 84)], c='#4daf4a', capsize=4)
plt.scatter(88.5,  np.median(ngc6300_beams_nn3), c='#a65628', label='3rd Nearest Neighbor')
plt.errorbar(88.5, [np.percentile(ngc6300_beams_nn3, 16),np.percentile(ngc6300_beams_nn3, 84)], c='#a65628', capsize=4)

plt.scatter(87.75,  np.median(ngc6300_match_beams_nn), c='#377eb8', marker='x', label='1st Nearest Neighbor')
plt.errorbar(87.75, [np.percentile(ngc6300_match_beams_nn, 16),np.percentile(ngc6300_match_beams_nn, 84)], c='#377eb8', capsize=4)
plt.scatter(88.2,  np.median(ngc6300_match_beams_nn2), c='#4daf4a', marker='x', label='2nd Nearest Neighbor')
plt.errorbar(88.2, [np.percentile(ngc6300_match_beams_nn2, 16),np.percentile(ngc6300_match_beams_nn2, 84)], c='#4daf4a', capsize=4)
plt.scatter(88.75,  np.median(ngc6300_match_beams_nn3), c='#a65628', marker='x', label='3rd Nearest Neighbor')
plt.errorbar(88.75, [np.percentile(ngc6300_match_beams_nn3, 16),np.percentile(ngc6300_match_beams_nn3, 84)], c='#a65628', capsize=4)

#Maxes
plt.scatter(90.6,  np.median(ngc3621_m_beams_nn), c='#377eb8', marker='s', label='1st Nearest Neighbor')
plt.errorbar(90.6, [np.percentile(ngc3621_m_beams_nn, 16),np.percentile(ngc3621_m_beams_nn, 84)], c='#377eb8', capsize=4)
plt.scatter(91.1,  np.median(ngc3621_m_beams_nn2), marker='s', c='#4daf4a', label='2nd Nearest Neighbor')
plt.errorbar(91.1, [np.percentile(ngc3621_m_beams_nn2, 16),np.percentile(ngc3621_m_beams_nn2, 84)], c='#4daf4a', capsize=4)
plt.scatter(91.6,  np.median(ngc3621_m_beams_nn3), c='#a65628', marker='s', label='3rd Nearest Neighbor')
plt.errorbar(91.6, [np.percentile(ngc3621_m_beams_nn3, 16),np.percentile(ngc3621_m_beams_nn3, 84)], c='#a65628', capsize=4)
plt.scatter(87.6,  np.median(ngc6300_m_beams_nn), c='#377eb8', marker='s', label='1st Nearest Neighbor')
plt.errorbar(87.6, [np.percentile(ngc6300_m_beams_nn, 16),np.percentile(ngc6300_m_beams_nn, 84)], c='#377eb8', capsize=4)
plt.scatter(88.1,  np.median(ngc6300_m_beams_nn2), c='#4daf4a', marker='s', label='2nd Nearest Neighbor')
plt.errorbar(88.1, [np.percentile(ngc6300_m_beams_nn2, 16),np.percentile(ngc6300_m_beams_nn2, 84)], c='#4daf4a', capsize=4)
plt.scatter(88.6,  np.median(ngc6300_m_beams_nn3), c='#a65628', marker='s', label='3rd Nearest Neighbor')
plt.errorbar(88.6, [np.percentile(ngc6300_m_beams_nn3, 16),np.percentile(ngc6300_m_beams_nn3, 84)], c='#a65628', capsize=4)

#AXIS & LABELS
plt.legend(handles=[first, second, third, centers, matched, maxes])
plt.xlabel('NGC 6300                  |                  NGC 3621', fontsize=15)
plt.ylabel('Distance to Nearest Neighbor (beams)')
#plt.ylim(0,np.max(merr[1,2]+5))
plt.xticks([])
plt.title('Max vs. Center', fontsize=16)
plt.savefig('max_v_center.pdf')
plt.show()
plt.close()
