from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import pandas as pd

cat60 = pd.read_csv('ngc3621_60pc_cloud_stats.csv')
cat90 = pd.read_csv('ngc3621_90pc_cloud_stats.csv')
cat120 = pd.read_csv('ngc3621_120pc_cloud_stats.csv')
cat150 = pd.read_csv('ngc3621_150pc_cloud_stats.csv')

first60 = cat60['beam_sep_nn']
second60 = cat60['beam_sep_nn2']
third60 = cat60['beam_sep_nn3']

first90 = cat90['beam_sep_nn']
second90 = cat90['beam_sep_nn2']
third90 = cat90['beam_sep_nn3']

first120 = cat120['beam_sep_nn']
second120 = cat120['beam_sep_nn2']
third120 = cat120['beam_sep_nn3']

first150 = cat150['beam_sep_nn']
second150 = cat150['beam_sep_nn2']
third150 = cat150['beam_sep_nn3']


plt.scatter(second60, third60, label='60 pc', alpha = 0.25)
plt.scatter(second90, third90, label='90 pc', alpha = 0.25)
plt.scatter(second120, third120, label='120 pc', alpha = 0.25)
plt.scatter(second150, third150, label='150 pc', alpha = 0.25)
plt.xlabel('Distance to Second Nearest Neighbor (beams)')
plt.ylabel('Distance to Third Nearest Neighbor (beams)')
plt.xlim(0,12)
plt.ylim(0,15)
plt.legend()
plt.savefig('second_vs_third_beam.pdf')
plt.close()
