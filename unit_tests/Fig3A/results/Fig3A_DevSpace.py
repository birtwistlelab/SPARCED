# %%
import os 
import sys
import jdata as jd
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

os.chdir('/home/jonah/Desktop/SPARCED/unit_tests/Fig3A/results')
data = jd.load('Fig3A.json')

# %%
ligInds = [155,156,157,158,159,160,161] # E, H, HGF, P, F, I, INS
ligConc = [0.001,0.01,0.1,1,10,100,1000]
ResInds = [162,166,167,170,171,172,173,174] # E1, E3, E4, Met, Pr, Fr, Ir, Isr
RecConc = 60.8754593

conds = {}

SPInds = np.arange(188,194) #EE1, HE3, HE4, HGF_Met, PPr, FFr, EE1E2
SPPInds = np.arange(194,214) #EE1E2, 
SPSPInds = np.arange(214,225) #EE1E2E3, EE1HE3, EE1HE4, HE3HE3, HE3HE4, HE4HE4, HGF_Met_HGF_Met, PPrPPr, FFrFFr, IIrIrI, INS_Isr_Isr_INS, FFE1_ppERK
for condition in data:
      ns_prod_sum = sum(data[condition]['xoutS'][-1, SPInds]) + sum(data[condition]['xoutS'][-1, SPPInds])\
            + 2.0*sum(data[condition]['xoutS'][-1, SPSPInds])
      conds[condition] = ns_prod_sum
            # print(f'iteration {ligand} {str(conc)} {ns_prod_sum}')      

egf_results = list(conds.values())[0:7]
h_results = list(conds.values())[7:14]
# hgf_results = conds[16:24]
hgf_results = list(conds.values())[14:21]
# p_results = conds[24:32]
p_results = list(conds.values())[21:28]
# f_results = conds[32:40]
f_results = list(conds.values())[28:35]
# i_results = conds[40:48]
i_results = list(conds.values())[35:42]
# ins_results = conds[48:56]
ins_results = list(conds.values())[42:49]

list_of_results = [egf_results, h_results, hgf_results, p_results, f_results, i_results, ins_results]
list_of_ligands = ['EGF', 'H', 'HGF', 'P', 'F', 'I', 'INS']

print(h_results)

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# Define the function you want to fit
def func(x, a, b, c):
    return a * np.power(x, b) / (np.power(c, b) + np.power(x, b))

# Calculate the number of rows based on the number of results
num_results = len(list_of_results)
num_columns = 3
num_rows = -(-num_results // num_columns)  # Ceiling division to ensure all results are shown

# Create a figure and a grid of subplots
fig, axes = plt.subplots(nrows=num_rows, ncols=num_columns, figsize=(10, 2.5 * num_rows))

# Loop through each lig_result and plot on separate subplots
for i, lig_result in enumerate(list_of_results):
    x_data = ligConc
    y_data = lig_result

    # Fit the curve
    params, covariance = curve_fit(func, x_data, y_data)
    a_fit, b_fit, c_fit = params
    y_fit = func(x_data, a_fit, b_fit, c_fit)

    # Plot on the i-th subplot
    row = i // num_columns
    col = i % num_columns
    ax = axes[row, col] if num_rows > 1 else axes[col]
    ax.scatter(x_data, y_data, label='Data with Noise')
    ax.plot(x_data, y_fit, label=f'Fitted Curve: a={a_fit:.2f}, b={b_fit:.2f}, c={c_fit:.2f}', color='red')
    ax.set_title(list_of_ligands[i])
    ax.tick_params(axis='y', labelsize=14)
    # ax.set_xlabel('x')
    # ax.set_ylabel('y')
    ax.set_xscale('log', base=10)
    # ax.legend()
# fig.xlabel('Ligand Concentration (nM)')
# Remove empty subplots if needed
for i in range(num_results, num_rows * num_columns):
    row = i // num_columns
    col = i % num_columns
    if num_rows > 1:
        fig.delaxes(axes[row, col])
    else:
        fig.delaxes(axes[col])

# Adjust layout to prevent clipping of titles and labels
plt.tight_layout()

# Show the entire figure
plt.show()

# %%



