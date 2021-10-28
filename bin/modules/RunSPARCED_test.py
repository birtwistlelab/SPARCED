import libsbml
import importlib
import amici
import numpy as np
import pandas as pd
import os


from modules.SGEmodule_test import SGEmodule_test
from modules.RunPrep import RunPrep

def RunSPARCED_test(flagD,th,spdata,Vn,Vc,model,wd,omics_input='OmicsData.txt',genereg_input='GeneReg.txt'):
    np.random.seed()
    ts = 30 # time-step to update mRNA numbers
    NSteps = int(th*3600/ts)
    tout_all = np.arange(0,th*3600+1,30)    
    mpc2nM_Vc = (1E9/(Vc*6.023E+23))
    
    genedata, mExp_mpc, GenePositionMatrix, AllGenesVec, kTCmaxs, kTCleak, kTCleak2, kGin_1, kGac_1, kTCd, TARs0, tcnas, tcnrs, tck50as, tck50rs, spIDs, mrna_idx = RunPrep(flagD,Vn,model,wd,omics_input,genereg_input)
    
    if len(spdata)==0:
        spdata0 = pd.read_csv(os.path.join(wd,'input_files','Species.txt'),header=0,index_col=0,sep="\t")
        spdata = np.float(spdata0.values[:,1])
    xoutS_all = np.zeros(shape=(NSteps+1,len(spdata)))
    xoutS_all[0,:] = spdata # 24hr time point     
    
    if len(genedata)==0:
        genedata = genedata0
    xoutG_all = np.zeros(shape=(NSteps+1,len(genedata)))
    xoutG_all[0,:] = genedata
    
    solver = model.getSolver() # Create solver instance
    solver.setMaxSteps = 1e10
    
    Nb_all = []
    Nd_all = []
    ac2in_all = []
    in2ac_all = []
    
    for qq in range(NSteps):
        genedata,xmN,AllGenesVec, sum_Nb, sum_Nd, sum_ac2in, sum_in2ac = SGEmodule_test(flagD,ts,xoutG_all[qq,:],xoutS_all[qq,:],Vn,Vc,kTCmaxs,kTCleak,kTCd,AllGenesVec,GenePositionMatrix,kGin_1,kGac_1,tcnas,tck50as,tcnrs,tck50rs,spIDs,mrna_idx)
        xoutS_all[qq,mrna_idx] = np.dot(xmN,mpc2nM_Vc)
        model.setInitialStates(xoutS_all[qq,:])
        rdata = amici.runAmiciSimulation(model, solver)  # Run simulation
        xoutS_all[qq+1,:] = rdata['x'][-1,:]
        xoutG_all[qq+1,:] = genedata
        
        Nb_all.append(sum_Nb)
        Nd_all.append(sum_Nd)
        ac2in_all.append(sum_ac2in)
        in2ac_all.append(sum_in2ac)
        
        
        if rdata['x'][-1,103] < rdata['x'][-1,105]:
            print('Apoptosis happened')
            break
    xoutS_all = xoutS_all[~np.all(xoutS_all == 0, axis=1)]
    xoutG_all = xoutG_all[~np.all(xoutG_all == 0, axis=1)]
    tout_all = tout_all[0:len(xoutS_all)]
    
    Nb_all = np.array(Nb_all)
    Nd_all = np.array(Nd_all)
    ac2in_all = np.array(ac2in_all)
    in2ac_all = np.array(in2ac_all)
    
    return xoutS_all, xoutG_all, tout_all, Nb_all, Nd_all, ac2in_all, in2ac_all