#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 19:21:31 2022

@author: arnab
"""
import plotly.figure_factory as ff
import numpy as np
import matplotlib.pyplot as plt
import plotly.io as pio
pio.renderers.default='browser'



np.random.seed(1)

X = np.random.rand(15,12)

fig = ff.create_dendrogram(X)

fig.update_layout(width = 800, height = 500)
fig.show()

#%%

import plotly.figure_factory as ff

import numpy as np

X = np.random.rand(10, 12)
names = ['Jack', 'Oxana', 'John', 'Chelsea', 'Mark', 'Alice', 'Charlie', 'Rob', 'Lisa', 'Lily']
fig = ff.create_dendrogram(X, orientation='left', labels=names)
fig.update_layout(width=800, height=800)
fig.show()


#%% newick tree sequence

# ((((C:10,G:10)E:22,(H:10,J:10)K:22)L:22,((C:10,G:10)E:22,(H:10,J:10)K:22)L:22)M:22,(((C:10,G:10)E:22,(H:10,J:10)K:22)L:22,((C:10,G:10)E:22,(H:10,J:10)K:22)L:22)M:22,(((C:10,G:10)E:22,(H:10,J:10)K:22)L:22,((C:10,G:10)E:22,(H:10,J:10)K:22)L:22)M:22,(((C:10,G:10)E:22,(H:10,J:10)K:22)L:22,((C:10,G:10)E:22,(H:10,J:10)K:22)L:22)M:22,(((C:10,G:10)E:22,(H:10,J:10)K:22)L:22,((C:10,G:10)E:22,(H:10,J:10)K:22)L:22)M:22,(((C:10,G:10)E:22,(H:10,J:10)K:22)L:22,((C:10,G:10)E:22,(H:10,J:10)K:22)L:22)M:22);