
import pandas as pd
import numpy as np
import re
import libsbml
import os
import sys
import importlib
import amici
import amici.plotting

import matplotlib.pyplot as plt
import matplotlib as mpl


mpl.rcParams['figure.dpi'] = 300

#%% load model and input files

# load model

sbml_file = "SPARCED.xml"
model_name= sbml_file[0:-4]
model_output_dir = model_name
sys.path.insert(0, os.path.abspath(model_output_dir))
model_module = importlib.import_module(model_name)
model = model_module.getModel()



sbml_reader = libsbml.SBMLReader()
sbml_doc = sbml_reader.readSBML(sbml_file)
sbml_model = sbml_doc.getModel()
params_model = []
[params_model.append(p.getId()) for p in sbml_model.getListOfParameters()]


ratelaw_sheet = np.array([np.array(line.strip().split("\t")) for line in open('Ratelaws.txt')])
ratelaw_data = np.array([line[1:] for line in ratelaw_sheet[1:]])





#%% read initial parameters

gene_params = pd.read_csv(os.path.join('input_files','OmicsData_extended.txt'), sep=',', index_col=0, header=0)
model_genes = gene_params.index





gene_params['kTLd'] = np.log(2)/gene_params['Protein_half_life_lit_h']/3600

gene_params['kTLd'][np.isnan(gene_params['kTLd'].values)] = np.log(2)/gene_params['Protein_half_life_Schwan_h'][np.isnan(gene_params['kTLd'].values)]/3600

gene_params['kTLd'][np.isnan(gene_params['kTLd'].values)] = 0

gene_params['kTL_nat_cells'] = gene_params['Exp Protein']*gene_params['kTLd']/gene_params['Exp RNA']

gene_params['kTL_nat_cells'][np.isnan(gene_params['kTL_nat_cells']) | np.isinf(gene_params['kTL_nat_cells'])] = 0

gene_params['kTL_nat'] = gene_params['kTL_nat_cells']

gene_params['kTL_nat'][np.where(gene_params['kTL_nat']==0)[0]] = gene_params['kTLnatLit_s'][np.where(gene_params['kTL_nat']==0)[0]]

gene_params['kTL_nat'][np.where(gene_params['kTL_nat']==0)[0]] = gene_params['kTLnatSchwan_s'][np.where(gene_params['kTL_nat']==0)[0]]

gene_params['kTL_nat'][np.array([list(model_genes).index(x) for x in ['CCND1', 'CCND2', 'CCND3']])] = gene_params['kTL_nat'][np.array([list(model_genes).index(x) for x in ['CCND1', 'CCND2', 'CCND3']])].values * 5


Step1_mrna = pd.read_csv(os.path.join('input_files','initializer','Initializer.csv'),sep=',',squeeze=True,usecols=['Step1_mrna','Step1_mrna_mpc'],index_col='Step1_mrna')
Step1_mrna = Step1_mrna[Step1_mrna.index.notnull()]



for gene_symbol in Step1_mrna.index:
    gene_params.loc[gene_symbol,'Exp RNA'] = Step1_mrna[gene_symbol]

#%% test - read reaction id from sbml

# reactions_all = [sbml_model.getListOfReactions().get(r).getId() for r in range(len(sbml_model.getListOfReactions()))]

reactions_all = [ratelaw_sheet[i][0] for i in range(len(ratelaw_sheet))]

vTL_pattern = re.compile("vTL\d+")

reactions_vTL = list(filter(vTL_pattern.match, reactions_all))




#%% parameter IDs

numberofgenes = len(model_genes)
S_PARCDL = pd.read_csv(os.path.join('input_files',"StoicMat.txt"), header=0, index_col=0, sep='\t')
# S_TL = S_PARCDL.iloc[:,2:(2+numberofgenes)]
S_TL = S_PARCDL.loc[:,S_PARCDL.columns.isin(reactions_vTL)]

ObsMat = pd.read_csv(os.path.join('input_files','Observables.txt'), header=0, index_col=0, sep='\t')
NumObs = len(ObsMat.columns)


kTL_id = []
kTL_default = []
k50E_id = []
k50E_values = []
kTL_genes = []


Step1_par = pd.Series()

for rowNum, ratelaw in enumerate(ratelaw_data):
    
    if "kTL" in str(ratelaw[1]):
        kTL_i = "k"+str(rowNum+1)+"_1"
        kTL_id.append(kTL_i)
        kTL_gene = re.search(r'm_\w+', ratelaw[1]).group()[2:]
        kTL_genes.append(kTL_gene)
        params_i = np.array(list(ratelaw[2:]))
        params_i = np.array([float(params_i[k]) for k in range(len(params_i))])
        params_i = params_i[~np.isnan(params_i)]
        kTL_value = params_i[0]
        kTL_default.append(float(kTL_value))
        a = len(params_i)
        if "EIF4E" in str(ratelaw[1]):
            
            if 'MDM2pro' in str(ratelaw[1]):
                k50E_id_i = "k"+str(rowNum+1)+"_3"
            elif 'CCND1' in str(ratelaw[1]) or 'CCND2' in str(ratelaw[1]) or 'CCND3' in str(ratelaw[1]):
                k50E_id_i = "k"+str(rowNum+1)+"_5"
            else:
                k50E_id_i = "k"+str(rowNum+1)+"_2"
            k50E_id.append(k50E_id_i)
            k50E_value_i = params_i[-1]
            k50E_values.append(k50E_value_i)   
    
    if ratelaw_sheet[rowNum][0] == 'vbR':
        kbR0_id = "k"+str(rowNum)+"_1"
        kbRi_id = "k"+str(rowNum)+"_2"
        
    if ratelaw_sheet[rowNum][0] == 'vdR':
        kdR0_id = "k"+str(rowNum)

    if ratelaw_sheet[rowNum][0] == 'vA77':
        kA77_id = "k"+str(rowNum)
        Step1_par = Step1_par.append(pd.Series(kA77_id,index=['kA77']))
    if ratelaw_sheet[rowNum][0] == 'vA87':
        kA87_id = "k"+str(rowNum)
        Step1_par = Step1_par.append(pd.Series(kA87_id,index=['kA87']))
    if ratelaw_sheet[rowNum][0] == 'vC104':
        kC82_id = "k"+str(rowNum)
        Step1_par = Step1_par.append(pd.Series(kC82_id,index=['kC82']))
    if ratelaw_sheet[rowNum][0] == 'vD66':
        kDDbasal_id = "k"+str(rowNum)+"_1"
        Step1_par = Step1_par.append(pd.Series(kDDbasal_id,index=['kDDbasal']))
        



#%%




obs2exclude = ObsMat.columns[:26]
# obs2include = ObsMat.columns[26:]

kTLest = gene_params['kTL_nat'].values




mExp_mpc = gene_params['Exp RNA'].copy()

cell_params = pd.read_csv(os.path.join('input_files',"Compartments.txt"), header=0, index_col=0, sep='\t')
Vc = cell_params.loc['Cytoplasm','volume']
Vn = cell_params.loc['Nucleus','volume']
Vm = cell_params.loc['Mitochondrion','volume']
Ve = cell_params.loc['Extracellular','volume']
volumeofcell = Vc + Vn



fileSpecies = 'Species.txt' # input
ICf = pd.read_csv(os.path.join('input_files',fileSpecies),header=0,index_col=0,sep='\t')


mrna_id = []
mrna_filter = filter(lambda a: a.startswith('m_'), list(ICf.index))
for m in mrna_filter:
    mrna_id.append(m)





VxPARCDL = ICf.loc[:,'compartment'].copy()


for i in range(len(VxPARCDL)):
    if VxPARCDL[i] == "Cytoplasm":
        VxPARCDL[i] = Vc
    elif VxPARCDL[i] == "Nucleus":
        VxPARCDL[i] = Vn
    elif VxPARCDL[i] == "Mitochondrion":
        VxPARCDL[i] = Vm
    elif VxPARCDL[i] == "Extracellular":
        VxPARCDL[i] = Ve


VxTL = np.ones(numberofgenes)*Vc
for i in range(np.shape(S_TL)[1]):
    if len(np.nonzero(S_TL.values[:,i])[0]) != 0:
        obs_ind = int(np.nonzero(S_TL.values[:,i])[0])
        VxTL[i] = VxPARCDL[obs_ind]
        

xp_mpc = gene_params['kTL_nat'].values*mExp_mpc.values/gene_params['kTLd'].values

xp_mpc[np.isnan(xp_mpc)] = 0

pExp_nM = xp_mpc*1e9/(VxTL*6.023e23)
x0PARCDL = np.matmul(S_TL.values,pExp_nM)
x0PARCDL = pd.Series(x0PARCDL)
x0PARCDL.index = S_TL.index
mpc2nmcf_Vc=1E9/(Vc*6.023E+23)
mExp_nM=mExp_mpc*mpc2nmcf_Vc



for m in mrna_id:
    x0PARCDL[m] = ICf.loc[m, 'IC_Xinitialized']


Step1sp = pd.read_csv(os.path.join('input_files','initializer','Initializer.csv'),sep=',',squeeze=True,usecols=['Step1_sp','Step1_val'],index_col='Step1_sp')

for sp in Step1sp.index:
    x0PARCDL[sp] = Step1sp[sp]



k50E_default = max(k50E_values)

kA77 = 0 #BIM*Bax
kA87 = 0 #C8 activation
kC82 = 0.06/3600 ##p21 degradation
kbRi = 0
kdR0 = 0
kbR0 = 0

#%% functions

x0 = x0PARCDL

def get_observables(xout, VxPARCDL, Vc):    
    #ObsMat = pd.read_excel("observables_mat_v4.csv", header=0, index_col=0)
    Obs = []  
    Vr = VxPARCDL/Vc
    for i in range(np.shape(ObsMat)[1]):
        Obs_i = np.sum(ObsMat.values[:,i]*xout*Vr.flatten())
        Obs.append(Obs_i)
    Obs = [0 if i <= 1e-6 else i for i in Obs]
    Obs = np.array(Obs)    
    return Obs

obs0 = get_observables(x0.values, VxPARCDL.values, Vc)
def array_editor(array, inds, x0):
    for k in range(len(inds)):
        array[inds[k]] = x0[k]


def fse(x,y):    
    fs_error = ((x-y)/x)    
    for i in range(len(fs_error)):
        if np.isnan(fs_error[i]):
            fs_error[i] = 0         
    return fs_error

def fe(x,y):    
    f_error = ((x-y)/x)    
    if np.isnan(f_error):
        f_error = 0
    if x == 0 and y == 0:
        f_error = 0         
    return f_error

def timecourse(species,rdata):    
    timeh = rdata['t']/3600
    species_ind = np.nonzero(S_PARCDL.index==species)[0][0]        
    x_t = rdata['x'][:,species_ind]
    plt.scatter(timeh,x_t)
    plt.ylabel(str(species))
    plt.xlabel('time(h)')
    plt.title('species timecourse')
    plt.show



def timecourse_obs(obs_name,rdata,obs0_def=obs0):
    
    timeh = rdata['t']/3600
    obs_ind = np.nonzero(ObsMat.columns==str(obs_name))[0][0]
    obs_t = rdata['y'][:,obs_ind]
    plt.figure()
    plt.scatter(timeh,obs_t)
    plt.axhline(y=obs0_def[obs_ind],color='r')
    plt.ylim(0,max(max(obs_t),obs0_def[obs_ind])*1.5)
    plt.ylabel(str(obs_name))
    plt.xlabel('time(h)')
    plt.title('obs timecourse')
    plt.show
    
    

#%% obs to exclude


kTL_mod = np.ones(numberofgenes)*0.25

def obs2gene (obs_name):
    gene = model_genes[np.nonzero(np.matmul(ObsMat.values[:,list(ObsMat.columns).index(obs_name)],S_TL.values))[0]]
    return gene

def obs2gene_i (obs_name):
    gene_i = np.nonzero(np.matmul(ObsMat.values[:,list(ObsMat.columns).index(obs_name)],S_TL.values))[0]
    return gene_i


for m in obs2exclude:
    kTL_mod[obs2gene_i(m)] = 1.0
    kTLest[obs2gene_i(m)]
    for k in range(len(obs2gene_i(m))):
        kTLest[obs2gene_i(m)[k]] = model.getFixedParameterById(np.array(kTL_id)[obs2gene_i(m)][k])
    




#%% prep optimizer

ts = 1000*3600
model.setTimepoints(np.linspace(0,ts,1000))
model.setFixedParameterById(kA77_id, kA77)
model.setFixedParameterById(kA87_id, kA87)
model.setFixedParameterById(kC82_id, kC82)
model.setFixedParameterById(kbR0_id, kbR0)
model.setFixedParameterById(kbRi_id, kbRi)
model.setFixedParameterById(kdR0_id, kdR0)
model.setFixedParameterById('k316_1', 0.0) #kDDbasal


[model.setFixedParameterById(k50E_id[k],0) for k in range(len(k50E_id))]


x0 = x0PARCDL
model.setInitialStates(x0.values)


solver = model.getSolver()
solver.setMaxSteps = 1e10

kC173 = 1.1111e-4

kTL10_12 = kC173*17/sum(mExp_nM[9:12])
model.setFixedParameterById('k12_1',kTL10_12)
model.setFixedParameterById('k13_1',kTL10_12)
model.setFixedParameterById('k14_1',kTL10_12)


kTLest[9:12] = kTL10_12

kTL_initial = kTLest





def kTLadjustwhile(model,solver,x0, obs0, kTL_id, kTLest, kTL_mod, k50E_id, k50E_values, ObsMat, S_TL, flagE, flagR):
    
    if flagE == 0:
        [model.setFixedParameterById(k50E_id[k],0) for k in range(len(k50E_id))]
    elif flagE == 1:
        [model.setFixedParameterById(k50E_id[k],k50E_values[k]) for k in range(len(k50E_id))]
    
    if flagR == 0:
        model.setFixedParameterById(kbR0_id,0)
        model.setFixedParameterById(kbRi_id,0)
        model.setFixedParameterById(kdR0_id,0)
    elif flagR == 1:
        kbRi = 0.04
        kdR0 = 2.22e-6
        nR = 5
        k50R = 4.86        
        ps6ki = x0['pS6K']
        f1 = (ps6ki**nR)/(k50R**nR+ps6ki**nR)
        Rt = x0['Ribosome']
        kbR0 = Rt*kdR0 - kbRi*f1
        
        model.setFixedParameterById(kbR0_id, kbR0)
        model.setFixedParameterById(kbRi_id, kbRi)
        model.setFixedParameterById(kdR0_id, kdR0)
    
    model.setInitialStates(x0.values)
    
    [model.setFixedParameterById(kTL_id[k],kTLest[k]) for k in range(len(kTL_id))]
    
    m = len(ObsMat.columns)
    margin = 0.01
    
    while m!=0:
        rdata = amici.runAmiciSimulation(model,solver)
        obs1 = rdata['y'][-1]
        error_fe = (obs0 - obs1)/obs0
        kTLf_obs = np.ones(len(obs0))
        for i in range(len(error_fe)):
            if ObsMat.columns[i] in obs2exclude:
                kTLf_obs[i] = 1
            elif obs0[i] == 0:
                kTLf_obs[i] = 0        
            elif error_fe[i] > margin and ~np.isinf(error_fe[i]):
                kTLf_obs[i] = 1/(1-error_fe[i])
            elif error_fe[i] < -1 * margin and ~np.isinf(error_fe[i]):
                kTLf_obs[i] = 1/(1-error_fe[i])
            elif error_fe[i] > -1 * margin and error_fe[i] < margin:
                kTLf_obs[i] = 1
        
        kTLf = []
        for i in range(len(S_TL.columns)):
            a = np.nonzero(np.array(S_TL.iloc[:,i]))[0]
            if len(a) != 0:
                sp_ind = a[0]
                obs_ind = np.nonzero(np.array(ObsMat.iloc[sp_ind,:]))[0][0]
                kTLf.append(kTLf_obs[obs_ind])
            else:
                kTLf.append(1)
        kTLf = pd.Series(kTLf)
        kTLf = kTLf.transform(lambda x: 1 if np.isinf(x) or np.isnan(x) else x)
        kTLf = np.array(kTLf)
        
        kTLest = kTLest*(1+(kTLf-1)*kTL_mod)
        [model.setFixedParameterById(kTL_id[k],kTLest[k]) for k in range(len(kTL_id))]
        
        obs_notmatched = ObsMat.columns[~((error_fe > -1*margin) & (error_fe < margin) | (obs0==0))]
        obs_notmatched = obs_notmatched[~np.isin(obs_notmatched,obs2exclude)]
        m = len(obs_notmatched)
        
    kTLnew = kTLest
    [model.setFixedParameterById(kTL_id[k], kTLnew[k]) for k in range(len(kTL_id))]
    rdata_new = amici.runAmiciSimulation(model,solver)
    x_new = rdata_new['x']
    x1 = pd.Series(data=rdata_new['x'][-1], index=ObsMat.index)
    x1[x1.values<1e-6] = 0.0
    flagA = 0
    
    apop_def = x0['PARP']*.5
        
    parp_all = x_new[:,list(ObsMat.index).index('PARP')]
    
    if sum(parp_all < apop_def) != 0:
        flagA = 1    
    
    return kTLnew, rdata_new, x1, flagA


#%% adjust kTLs

kTLnew1, rdata_new, x1, flagA = kTLadjustwhile(model,solver,x0, obs0, kTL_id, kTLest, kTL_mod, k50E_id, k50E_values, ObsMat, S_TL, flagE=0, flagR=0)

kTLnew2, rdata_new, x2, flagA = kTLadjustwhile(model,solver,x1, obs0, kTL_id, kTLnew1, kTL_mod, k50E_id, k50E_values, ObsMat, S_TL, flagE=1, flagR=0)

kTLnew3, rdata_new, x3, flagA = kTLadjustwhile(model,solver,x2, obs0, kTL_id, kTLnew2, kTL_mod, k50E_id, k50E_values, ObsMat, S_TL, flagE=1, flagR=1)




#%% adjust Cd, p21

totalcyclinDfromdata = sum(pExp_nM[9:12])
totalp21fromdata = pExp_nM[25]

kC173 = 1.1111e-4
kC82 = 1.6667e-5*17

x_in = rdata_new['x'][-1]

th=0.001

ratio_cd=0.5
ratio_p21 = 0.5
model.setInitialStates(x_in)



kTL10_12 = kC173*17/sum(mExp_nM[9:12])
model.setFixedParameterById('k12_1',kTL10_12)
model.setFixedParameterById('k13_1',kTL10_12)
model.setFixedParameterById('k14_1',kTL10_12)


model.setFixedParameterById(kC82_id,kC82)

ratio_cd_list = []
kC173_list = []

cd_sp = np.argwhere(ObsMat.loc[:,'Cd'].values>0).flatten()
p21_sp = np.argwhere(ObsMat.loc[:,'p21'].values>0).flatten()


while (ratio_cd < (1-th) or ratio_cd > (1+th)) or (ratio_p21 < (1-th) or ratio_p21 > (1+th)):

    rdata_loop = amici.runAmiciSimulation(model,solver)
    
    ratio_cd = totalcyclinDfromdata/sum(rdata_loop['x'][-1][cd_sp])
    
    if ratio_cd < (1-th) or ratio_cd > (1+th):        
        f = 1 + (ratio_cd-1)*0.5
        kC173 = kC173*f
        kTL10_12 = kC173*17/sum(mExp_nM[9:12])
        model.setFixedParameterById('k12_1',kTL10_12)
        model.setFixedParameterById('k13_1',kTL10_12)
        model.setFixedParameterById('k14_1',kTL10_12)

    kC173_list.append(kC173)
    ratio_cd_list.append(ratio_cd)
    
    ratio_p21 = totalp21fromdata/sum(rdata_loop['x'][-1][p21_sp])
    
    if ratio_p21 < (1-th) or ratio_p21 > (1+th):  
        p = 1 + (ratio_p21-1)*0.5
        kC82 = kC82/p
        model.setFixedParameterById(kC82_id,kC82)
 
    
x4 = pd.Series(data=rdata_loop['x'][-1], index=ObsMat.index)
x4[x4.values<1e-6] = 0.0

kTLnew3[9:12] = kTL10_12
#%% adjust c8

kA77 = 3.162075e-9



model.setFixedParameterById(kA77_id, kA77)


kA87s = [1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]

for k in range(len(kA87s)):
    
    model.setFixedParameterById(kA87_id, kA87s[k])
    
    kTLnew4, rdata_loop, x5, flagA = kTLadjustwhile(model,solver,x4, obs0, kTL_id, kTLnew3, kTL_mod, k50E_id, k50E_values, ObsMat, S_TL, flagE=1, flagR=1)
    
    if flagA == 0:
        x5last = x5
        kTLnew4last = kTLnew4
    
    if flagA == 1:
        kA87 = kA87s[k-1]
        model.setFixedParameterById(kA87_id, kA87)
        break

x5 = x5last
kTLnew4 = kTLnew4last


x5[x5.values<1e-6] = 0.0

model.setInitialStates(x5.values)
[model.setFixedParameterById(kTL_id[k], kTLnew4[k]) for k in range(len(kTL_id))]



#%% adjust basal dna damage

BRCA2 = x5['BRCA2']
MSH6 = x5['MSH6']
MGMT = x5['MGMT']
Me = x5['Me']
Ma = x5['Ma']

fixdsb1 = model.getFixedParameterById('k313_1')
fixmsh = model.getFixedParameterById('k314_1')
fixmgmt = model.getFixedParameterById('k315_1')
kDDE = model.getFixedParameterById('k316_2')
kDEtop = model.getFixedParameterById('k316_4')
Etop = model.getFixedParameterById('k316_3')
kDnSP = model.getFixedParameterById('k316_5')
kDkmSP = model.getFixedParameterById('k316_6')
kDkmSS = model.getFixedParameterById('k267_3')
kDkmDS = model.getFixedParameterById('k264_3')

kDDbasal = 1e-6 # set this value manually

damageDSB_cycling = kDDbasal/(fixdsb1*BRCA2)
damageSSB_cycling = kDDbasal/(fixmsh*MSH6+fixmgmt*MGMT)

#%%
if damageDSB_cycling > kDkmDS:
    np.disp('ERROR --- DSB damage is too high, must reduce kDDbasal or increase strength of repair')
    
if damageSSB_cycling > kDkmSS:
    np.disp('ERROR --- SSB damage is too high, must reduce kDDbasal or increase strentgh of repair')


#%%

vdamage_on = (kDDbasal + kDDE*(Etop/(Etop+kDEtop)))*(((Me+Ma)**kDnSP)/(((Me+Ma)**kDnSP)+(kDkmSP**kDnSP)))

damageDSB = vdamage_on/(fixdsb1*BRCA2)
damageSSB = vdamage_on/(fixmsh*MSH6+fixmgmt*MGMT)

model.setFixedParameterById('k316_1', kDDbasal)

    

#%%

# Run simulation and get new steady state

x5['damageDSB'] = damageDSB
x5['damageSSB'] = damageSSB


rdata_new = amici.runAmiciSimulation(model,solver)

x6 = rdata_new['x'][-1]
x6[x6<1e-6] = 0
x6 = pd.Series(data=x6, index=ObsMat.index)



#%%

pcFos_cJun = x6['pcFos_cJun']
cMyc = x6['cMyc']
p53ac = x6['p53ac']
FOXOnuc = x6['FOXOnuc']
ppERKnuc = x6['ppERKnuc']
pRSKnuc = x6['pRSKnuc']
bCATENINnuc = x6['bCATENINnuc']

genereg = pd.read_csv(os.path.join('input_files','GeneReg.txt'), sep='\t', header=0, index_col=0)

numberofgenes = int(len(genereg.index))
numberofTARs = int(len(genereg.columns))

TAs = np.zeros([numberofgenes,numberofTARs])


#%%
mExp_mpc[5:9] = 17.0
mExp_mpc[12:25] = 17.0






kGsRead = pd.read_csv(os.path.join('input_files','OmicsData.txt'),header=0,index_col=0,sep='\t')
gExp_mpc = np.float64(kGsRead.loc[:,'Exp GCN'].values)
kGin = np.float64(kGsRead.loc[:,'kGin'].values)
kGac = np.float64(kGsRead.loc[:,'kGac'].values)
kTCmaxs = np.float64(kGsRead.loc[:,'kTCmaxs'].values)
kTCd = np.float64(kGsRead.loc[:,'kTCd'].values)

xgac_mpc_D = (kGac*gExp_mpc)/(kGin+kGac)

TARsRead = pd.read_csv('GeneReg.txt',header=0,index_col=0,sep="\t")
TARs0 = (TARsRead.values)

tcnas = np.ones((numberofgenes, numberofTARs))
tck50as = np.zeros((numberofgenes, numberofTARs))
tcnrs = np.ones((numberofgenes, numberofTARs))
tck50rs = np.zeros((numberofgenes, numberofTARs))

for qq in range(numberofgenes):
    for ww in range(numberofTARs):
        pars = TARs0[qq,ww].find(';')
        if pars>0:
            nH = np.float(TARs0[qq,ww][0:pars])
            kH = np.float(TARs0[qq,ww][pars+2::])
            if nH>0:
                tcnas[qq,ww] = nH
                tck50as[qq,ww] = kH
            else:
                tcnrs[qq,ww] = abs(nH)
                tck50rs[qq,ww] = kH

mpc2nmcf_Vn = 1.0E9/(Vn*6.023E+23)
# Convert to molecules per cell
tck50as = tck50as*(1/mpc2nmcf_Vn)
tck50rs = tck50rs*(1/mpc2nmcf_Vn)

spnames = [ele for ele in model.getStateIds()]
spIDs = []
for qq in range(numberofTARs):
    sps = spnames.index(TARsRead.columns[qq]) 
    spIDs.append(sps)

    
TARarr = np.array(x6[spIDs])
TAs = np.zeros((numberofgenes,numberofTARs))
TRs = np.zeros((numberofgenes,numberofTARs))
for qq in range(numberofTARs):
    TAs[tck50as[:,qq] > 0, qq] = TARarr[qq]
    TRs[tck50rs[:,qq] > 0, qq] = TARarr[qq]
TAs = TAs*(1.0/mpc2nmcf_Vn) # convert to mpc from nM
TAs.flatten()
TRs = TRs*(1.0/mpc2nmcf_Vn)
TRs.flatten() 

# make hills
aa = np.divide(TAs,tck50as)
TFa = np.power(aa,tcnas)
TFa[np.isnan(TFa)] = 0.0
bb = np.divide(TRs,tck50rs)
TFr = np.power(bb,tcnrs)
TFr[np.isnan(TFr)] = 0.0
hills = np.sum(TFa,axis=1)/(1 + np.sum(TFa,axis=1) + np.sum(TFr,axis=1))
# With AP1*cMYC exception:
hills[9:12] = np.multiply((TFa[9:12,0]/(1+TFa[9:12,0])),(TFa[9:12,1]/(1+TFa[9:12,1])))

# vTCd
vTCd= np.transpose(np.multiply(kTCd,mExp_mpc));
vTCd = np.squeeze(np.asarray(vTCd))

induced = np.multiply(np.multiply(xgac_mpc_D,kTCmaxs),hills)
induced = induced.flatten()

negativecheck = np.array(mExp_mpc,dtype=bool).astype(int)*(vTCd - induced)

i2c = np.nonzero(negativecheck<0)[0]

if len(i2c)!=0:
    np.disp('WARNING -- Some induction term exceed degradation terms')
    
leak = vTCd-induced
kTCleak_new=leak/xgac_mpc_D

kTCleak_new[np.isnan(kTCleak_new)] = 0
kTCleak_new[np.isinf(kTCleak_new)] = 0

#%%

