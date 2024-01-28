import os 
import sys 
import pandas as pd
import matplotlib.pyplot as plt
import jdata as jd

sys.path.append('/home/jonah/Desktop/SPARCED/unit_tests/src/')

os.chdir('/home/jonah/Desktop/SPARCED/unit_tests/Fig3E/petab_files')
data = jd.load('../results/Fig3E.json')

for condition in data:
    plt.plot(data[condition]['toutS']/3600, data[condition]['xoutS'][:, 2], label=condition)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xlabel('Time (hours)', fontsize=14)
    plt.ylabel('p53 active (mM)', fontsize=14)
    plt.title('Fig3E SPARCED_ERM Fig3E Unit Test Replication', fontsize=14)


plt.show()
plt.savefig('../results/Fig3E.png')
