#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import sys 
import jdata as jd

# os.chdir('/home/jrhuggi/projects/foreman/SPARCED/unit_tests/Fig4A/results/')

data = jd.load('Fig4A.json')


# In[2]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# In[3]:


def collect_the_dead(data):
    unique_conditions = set([stim for cells in data for stim in data[cells]])
    perished_cells = {stim: [] for stim in unique_conditions}
    for cell in data:
        for stim in data[cell]:
            perished_cell = np.argwhere(data[cell][stim]['xoutS'][:, 105]>100.0)
            if len(perished_cell) > 0: 
                death_point = perished_cell[0] #Instance where we rule apoptosis is irreversible
                time_of_death = data[cell][stim]['toutS'][death_point] # Timepoint to match death point. 
                perished_cells[stim].extend(time_of_death/3600)
            
    # Our returned array now has every cell listed, and the timepoints for when it died for each condition
    # Ex: condition1: cell 0: 48.5hours
    return perished_cells


# In[4]:


dict1 = collect_the_dead(data)


# In[5]:


merged_dict = {}

# for condition, cells in reorganized_dict.items():
for condition in dict1.keys():
    # Extract the base condition name (e.g., 'E_0_INS_0_k316_3_24')
    base_condition = condition.split('_24')[0]
    # print(base_condition)
    # If the base condition is not in merged_dict, create an empty list
    if base_condition not in merged_dict:
        merged_dict[base_condition] = []
    # Append the cells from the current condition to the base condition
    merged_dict[base_condition].extend(dict1[condition])


# In[6]:


def time_of_death(data, timepoints_of_interest):
    tod_timepoints = {}
    for condition in data:
        tod_timepoints[condition] = {}
        for timepoint in timepoints_of_interest:
            tod_timepoints[condition][timepoint] = []
            for tod in data[condition]:
                if tod < timepoint:
                    tod_timepoints[condition][timepoint].append(tod)
    return tod_timepoints


# In[11]:


dict2 = time_of_death(merged_dict, [24, 48, 72])


# In[12]:


def death_ratios(data, population_size):
    death_ratios = {}
    for condition in data:
        death_ratios[condition] = {}
        for timepoint in data[condition]:
            conditional_death_ratio = len(data[condition][timepoint]) / population_size
            death_ratios[condition][timepoint] = conditional_death_ratio
    # return death_ratios
    return death_ratios
death_ratio = death_ratios(dict2, 80)
death_ratio


# In[32]:


def plot_death_ratios(death_ratios):
    # Plot for each condition
    for condition, ratios in death_ratios.items():
        plt.figure(figsize=(3, 6))  # Create a new figure for each condition
        plt.bar(list(ratios.keys()), list(ratios.values()), width=20)
        plt.title(f'Condition: {condition}')
        plt.xlabel('Time Points')
        plt.ylabel('Death %(Scale 0-1)')
        plt.grid(True)
        plt.ylim(0,1)
        plt.xticks(list(ratios.keys()))

        # Show the plot
        plt.savefig(f'{condition}.png')

plot_death_ratios(death_ratio)

