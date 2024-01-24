# %%
import os
import sys

os.chdir('/home/jonah/Desktop/SPARCED/unit_tests/Fig2B/results')
sys.path.append('/home/jonah/Desktop/SPARCED/bin')
sys.path.append('/home/jonah/Desktop/SPARCED/unit_tests/src')

import numpy as np
import matplotlib.pyplot as plt
# from modules.runSPARCED import runSPARCED
import jdata as jd

data = jd.load('Fig2B.json')
data

# %%
plt.figure(figsize=(7, 4))
plt.plot(data['average_sim']['toutS']/3600.0,data['average_sim']['xoutG'][:,115],'b',linewidth=4)
# plt.xlabel('Time (hrs)')
# plt.ylabel('MAPK1 act genes', multialignment='center')
plt.grid(True)
plt.xlim(0, 24)
# plt.xticks(np.arange(0,24,step=4))
plt.xticks([])
plt.yticks(fontsize=16, weight='bold')
plt.savefig('Fig2B_g1.png')

plt.figure(figsize=(7, 4))
plt.plot(data['average_sim']['toutS']/3600.0,data['average_sim']['xoutG'][:,116],'r',linewidth=4)
# plt.xlabel('Time (hrs)')
# plt.ylabel('MAPK3 act genes', multialignment='center')
plt.grid(True)
plt.xlim(0, 24)
# plt.xticks(np.arange(0,24,step=4))
plt.xticks([])
plt.yticks(fontsize=16, weight='bold')
plt.savefig('Fig2B_g2.png')

# %%


plt.figure(figsize=(7, 4))
plt.plot(data['average_sim']['toutS']/3600.0,data['average_sim']['xoutS'][:,888],'b',linewidth=4)
# plt.xlabel('Time (hrs)')
# plt.ylabel('MAPK1 mRNAs', multialignment='center')
plt.grid(True)
plt.xlim(0, 24)
# plt.xticks(np.arange(0,24,step=4))
plt.xticks([])
plt.yticks(fontsize=16, weight='bold')
ax = plt.gca()
ax.yaxis.get_offset_text().set_fontsize(24)
ax.yaxis.get_offset_text().set_weight('bold')
plt.ticklabel_format(style='sci', scilimits=(-1, 0))
# plt.savefig('Fig2B_m1.png')

# plt.figure(figsize=(7, 4))
plt.plot(data['average_sim']['toutS']/3600.0,data['average_sim']['xoutS'][:,889],'r',linewidth=4)
# plt.xlabel('Time (hrs)')
# plt.ylabel('MAPK3 mRNAs', multialignment='center')
# plt.grid(True)
# plt.xlim(0, 24)
# plt.xticks(np.arange(0,24,step=4))
# plt.xticks([])
# plt.yticks(fontsize=14)
plt.savefig('Fig2B_m1_m2.png')



# %%


tERKinds = [150,154,225,226,227,674,675,715,716,717,728,729,732,735,736,737,739,741,749,755,756,757]

plt.figure(figsize=(7, 4))
plt.plot(data['average_sim']['toutS']/3600.0,data['average_sim']['xoutS'][:,tERKinds].sum(axis=1),'k',linewidth=4, color = 'purple')
# plt.xlabel('Time (hrs)')
# plt.ylabel('total ERK (nM)', multialignment='center')
plt.grid(True)
plt.xlim(0, 24)
plt.xticks(np.arange(0,25,step=4), fontsize=16, weight='bold')
plt.yticks(fontsize=16, weight='bold')
# plt.xticks([])
plt.savefig('Fig2B_tERK.png')

# %%



