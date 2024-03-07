import os 
import sys
import jdata as jd
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

os.chdir('/home/jonah/Desktop/SPARCED/unit_tests/Fig3A/results')
data = jd.load('Fig3A_alt.json')

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
      ns_prod_sum = sum(data[condition]['cell 0']['xoutS'][-1, SPInds]) + sum(data[condition]['cell 0']['xoutS'][-1, SPPInds])\
            + 2.0*sum(data[condition]['cell 0']['xoutS'][-1, SPSPInds])
      conds[condition] = ns_prod_sum
            # print(f'iteration {ligand} {str(conc)} {ns_prod_sum}')      

egf_results = list(conds.values())[0:7]
print(egf_results)
h_results = list(conds.values())[7:14]
print(h_results)
hgf_results = list(conds.values())[14:21]
print(hgf_results)
p_results = list(conds.values())[21:28]
print(p_results)
f_results = list(conds.values())[28:35]
print(f_results)
i_results = list(conds.values())[35:42]
print(i_results)
ins_results = list(conds.values())[42:49]
print(ins_results)

list_of_results = [egf_results, h_results, hgf_results, p_results, f_results, i_results, ins_results]
list_of_ligands = ['EGF', 'H', 'HGF', 'P', 'F', 'I', 'INS']

print(h_results)

from scipy.optimize import curve_fit

# Define the function you want to fit
def func(x, a, b, c):
    return a * np.power(x, b) / (np.power(c, b) + np.power(x, b))

# Calculate the number of rows based on the number of results
num_results = len(list_of_results)
num_columns = 3
num_rows = -(-num_results // num_columns)  # Ceiling division to ensure all results are shown


# Set global plot attributes
plt.rc('font', size=16, weight='bold')  # Set fontsize and font weight
plt.rc('axes', titlesize=16)             # Set title fontsize
# plt.rc('axes', labelsize=12)             # Set label fontsize
plt.rc('xtick', labelsize=16)            # Set x-axis tick label fontsize
plt.rc('ytick', labelsize=16)            # Set y-axis tick label fontsize
plt.rc('legend', fontsize=16)
# Create a figure and a grid of subplots
fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(10, 7.5))

hill_curve_data = {}

# Loop through each lig_result and plot on separate subplots
for i, lig_result in enumerate(list_of_results):
    hill_curve_data[list_of_ligands[i]] = {}
    x_data = ligConc
    y_data = lig_result
    # Fit the curve
    params, covariance = curve_fit(func, x_data, y_data)
    a_fit, b_fit, c_fit = params
    y_fit = func(x_data, a_fit, b_fit, c_fit)

    hill_curve_data[list_of_ligands[i]]['x'] = x_data
    hill_curve_data[list_of_ligands[i]]['y'] = y_data
    hill_curve_data[list_of_ligands[i]]['fit'] = y_fit
    hill_curve_data[list_of_ligands[i]]['hill coefficient'] = b_fit


axes[0, 0].scatter(hill_curve_data['EGF']['x'], hill_curve_data['EGF']['y'])
axes[0, 0].plot(hill_curve_data['EGF']['x'], hill_curve_data['EGF']['fit'], \
    label=f'n={hill_curve_data["EGF"]["hill coefficient"]:.2f}', color='blue')
axes[0, 0].set_title('EGF')
axes[0, 0].set_xscale('log', base=10)
axes[0, 0].legend(frameon=False)

axes[0, 1].scatter(hill_curve_data['H']['x'], hill_curve_data['H']['y'], color='black')
axes[0, 1].plot(hill_curve_data['H']['x'], hill_curve_data['H']['fit'], \
    label=f'n={hill_curve_data["H"]["hill coefficient"]:.2f}', color='black')
axes[0, 1].set_title('H')
axes[0, 1].set_xscale('log', base=10)
axes[0, 1].legend(frameon=False)

axes[0, 2].scatter(hill_curve_data['HGF']['x'], hill_curve_data['HGF']['y'], color='black')
axes[0, 2].plot(hill_curve_data['HGF']['x'], hill_curve_data['HGF']['fit'], \
    label=f'n={hill_curve_data["HGF"]["hill coefficient"]:.2f}', color='black')
axes[0, 2].set_title('HGF')
axes[0, 2].set_xscale('log', base=10)
axes[0, 2].legend(frameon=False)

axes[1, 0].scatter(hill_curve_data['P']['x'], hill_curve_data['P']['y'],color='red')
axes[1, 0].plot(hill_curve_data['P']['x'], hill_curve_data['P']['fit'], \
    label=f'n={hill_curve_data["P"]["hill coefficient"]:.2f}', color='red')
axes[1, 0].set_title('P')
axes[1, 0].set_xscale('log', base=10)
axes[1, 0].legend(frameon=False)

axes[1, 1].scatter(hill_curve_data['F']['x'], hill_curve_data['F']['y'])
axes[1, 1].plot(hill_curve_data['F']['x'], hill_curve_data['F']['fit'], \
    label=f'n={hill_curve_data["F"]["hill coefficient"]:.2f}', color='blue')
axes[1, 1].set_title('F')
axes[1, 1].set_xscale('log', base=10)
axes[1, 1].legend(frameon=False)

axes[1, 2].scatter(hill_curve_data['I']['x'], hill_curve_data['I']['y'])
axes[1, 2].plot(hill_curve_data['I']['x'], hill_curve_data['I']['fit'], \
    label=f'n={hill_curve_data["I"]["hill coefficient"]:.2f}', color='blue')
axes[1, 2].set_title('I')
axes[1, 2].set_xscale('log', base=10)
axes[1, 2].legend(frameon=False)

axes[2, 0].scatter(hill_curve_data['INS']['x'], hill_curve_data['INS']['y'])
axes[2, 0].plot(hill_curve_data['INS']['x'], hill_curve_data['INS']['fit'], \
    label=f'n={hill_curve_data["INS"]["hill coefficient"]:.2f}', color='blue')
axes[2, 0].set_title('INS')
axes[2, 0].set_xscale('log', base=10)
axes[2, 0].legend(frameon=False)


fig.delaxes(axes[2, 1])
fig.delaxes(axes[2, 2])
plt.tight_layout()

# Show the entire figure
plt.savefig('Fig3A.png')


