from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import pandas as pd

cat60 = pd.read_csv('ngc3621_60pc_cloud_stats.csv')
cat90 = pd.read_csv('ngc3621_90pc_cloud_stats.csv')
cat120 = pd.read_csv('ngc3621_120pc_cloud_stats.csv')
cat150 = pd.read_csv('ngc3621_150pc_cloud_stats.csv')

first60 = cat60['min_dist']
second60 = cat60['min_dist2nd']
third60 = cat60['min_dist3rd']

first90 = cat90['min_dist']
second90 = cat90['min_dist2nd']
third90 = cat90['min_dist3rd']

first120 = cat120['min_dist']
second120 = cat120['min_dist2nd']
third120 = cat120['min_dist3rd']

first150 = cat150['min_dist']
second150 = cat150['min_dist2nd']
third150 = cat150['min_dist3rd']


plt.scatter(second60, third60, label='60 pc', alpha = 0.25)
plt.scatter(second90, third90, label='90 pc', alpha = 0.25)
plt.scatter(second120, third120, label='120 pc', alpha = 0.25)
plt.scatter(second150, third150, label='150 pc', alpha = 0.25)
plt.xlabel('Distance to Second Nearest Neighbor (pc)')
plt.ylabel('Distance to Third Nearest Neighbor (pc)')
plt.xlim(0,1100)
plt.ylim(0,1200)
plt.legend()
plt.savefig('second_vs_third_pc.pdf')
plt.close()
