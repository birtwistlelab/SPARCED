import libsbml
import importlib
import amici
import numpy as np
import pandas as pd
import os

from modules.SGEmodule import SGEmodule
from modules.RunPrep import RunPrep

def RunSPARCED(flagD,th,spdata,genedata,sbml_file,model):
    ts = 30 # time-step to update mRNA numbers
    NSteps = int(th*3600/ts)
    tout_all = np.arange(0,th*3600+1,ts) 
    
    cd = os.getcwd()
    wd = os.path.dirname(cd)
    input_path = os.path.join(wd,'input_files')
    
    # Read-in the model SBML to get compartmental volumes (used to convert nM to mpc and vice versa)
    sbml_reader = libsbml.SBMLReader()
    sbml_doc = sbml_reader.readSBML(os.path.join(wd,sbml_file))
    sbml_model = sbml_doc.getModel()
    Vc = sbml_model.getCompartment(0).getVolume() # Should be the index for Cytoplasm
    Vn = sbml_model.getCompartment(2).getVolume() # Should be the index for Nuclues
    mpc2nM_Vc = (1E9/(Vc*6.023E+23))
    splist = list(model.getStateIds())
    if len(spdata)==0: # if no initial condition values are supplied, use the input file information
        spdata0 = pd.read_csv(os.path.join(input_path,'Species.txt'),header=0,index_col=0,sep="\t")
        spdata = np.float(spdata0.values[:,1])  
    
    # calculate 
    genedata, GenePositionMatrix, AllGenesVec, kTCmaxs, kTCleak, kGin_1, kGac_1, kTCd, TARs0, tcnas, tcnrs, tck50as, tck50rs, spIDs = RunPrep(flagD,Vn,model)
    
    
    xoutS_all = np.zeros(shape=(NSteps+1,len(spdata)))
    
    # xoutS_all = np.zeros(shape=(1,len(splist)))
    xoutS_all[0,:] = spdata # 24hr time point
    # xoutG_all = np.zeros(shape=(1,len(genedata)))
    xoutG_all = np.zeros(shape=(NSteps+1,len(genedata)))
    xoutG_all[0,:] = genedata  
    
    solver = model.getSolver() # Create solver instance
    solver.setMaxSteps = 1e10
    
    mRNAIndDs = [ind for ind, ele in enumerate(splist) if 'm_' in ele] # find the indeces for mRNA species
    mRNAIndDs = mRNAIndDs[1:]
    PARPind = [ind for ind,ele in enumerate(splist) if ele in {'PARP'}] # find the index for PARP
    cPARPind = [ind for ind,ele in enumerate(splist) if ele in {'cPARP'}] # find the index for cleaved-PARP (used to decide for apoptosis) 
    n_sp = len(splist)
    # Run 30sec (ts) simulations until final th is reached:
    for qq in range(NSteps): 
        # Call the function (based on the current state of the model species) for gene in/activation and mRNA birth/death events.   
        # Stochastic sampling if the flagD==0, deterministic calculations if flagD==1:
        genedata,xmN,AllGenesVec = SGEmodule(flagD,ts,xoutG_all[qq,:],xoutS_all[qq,:],Vn,Vc,kTCmaxs,kTCleak,kTCd,AllGenesVec,GenePositionMatrix,kGin_1,kGac_1,tcnas,tck50as,tcnrs,tck50rs,spIDs,mRNAIndDs[0])
        # mRNA species values are updated every 30sec, for the next 30sec simulation:
        xoutS_all[qq,mRNAIndDs] = np.dot(xmN,mpc2nM_Vc) 
        # set the new ICs:
        model.setInitialStates(xoutS_all[qq,:]) 
        # Run the simulation:
        rdata = amici.runAmiciSimulation(model, solver)  
        # Store the end point as next 30sec time-point:
            
        xoutS_all[qq+1,:] = rdata._swigptr.x[-n_sp:]
        # xoutS_all = np.vstack([xoutS_all, rdata['x'][-1,:]]) 
        rdata = None
        # Store active/inactive gene states:
        xoutG_all[qq+1,:] = genedata
        # xoutG_all = np.vstack([xoutG_all, genedata]) 
        # check for cell death:
        if xoutS_all[-1,PARPind] < xoutS_all[-1,cPARPind]: 
            print('Apoptosis happened')
            break
    # Finalize the species concentration trajectories (output at every 30sec):
    xoutS_all = xoutS_all[~np.all(xoutS_all == 0, axis=1)] 
    # Finalize the gene state trajectories (output at every 30sec):
    xoutG_all = xoutG_all[~np.all(xoutG_all == 0, axis=1)] 
    # The time points (in seconds):
    tout_all = tout_all[0:len(xoutS_all)]
    
    return xoutS_all, xoutG_all, tout_all
