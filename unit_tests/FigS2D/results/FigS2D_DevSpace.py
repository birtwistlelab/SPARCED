import os
import sys
import jdata as jd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
os.chdir('/home/jonah/Desktop/SPARCED/unit_tests/FigS2D/scripts')

sys.path.insert(0, '/home/jonah/Desktop/SPARCED/unit_tests/src')

data = jd.load('../results/FigS2D.json')

Vc = 5.25E-12

mpc2nmcf_Vc = 1.0E9/(Vc*6.023E+23)

plt.figure(figsize=(7, 4))
yy = np.array(data['Ribosome_mpc']['Ribosome_total']['xoutS']*(1/mpc2nmcf_Vc))
tt = (data['Ribosome_mpc']['Ribosome_total']['toutS']/3600)
plt.plot(tt, yy,linewidth=4, color='black')
plt.xlim([0, 24])
plt.ylabel('Ribosome (mpc)')
plt.xlabel('Time (hours)')
plt.savefig('../results/FigS2D_v1.png')
