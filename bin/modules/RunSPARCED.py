import libsbml
import importlib
import amici
import numpy as np

from modules.SGEmodule import SGEmodule
from modules.RunPrep import RunPrep

def RunSPARCED(flagD,th,species_initializations,genedata,Vn,Vc,model):
    # Inputs:
    # flagD = deterministic (1) or stochastic (0) simulation
    # ts = time step
    # genedata = 1D array of latest gene-expression module concentrations
    # spdata = 1D array of latest species concentrations
    # Vn = nuclear volume
    # Vc = cytoplasmic volume
    # model = The SBML model instance to run

    # Outputs:
    # xoutS_all = Concentrations of all species
    # xoutG_all = Numbers of active and inactive genes
    # tout_all = Time points of the outputs

    ts = 30 # time-step to update mRNA numbers
    NSteps = int(th*3600/ts)
    tout_all = np.arange(0,th*3600+1,30)
    mpc2nM_Vc = (1E9/(Vc*6.023E+23))

    genedata0,mRNA_mpc,GenePositionMatrix,AllGenesVec,kTCmaxs,kTCleak,kTCleak2,kGin_1,kGac_1,kTCd,TAs0,TRs0,tcnas,tcnrs, \
    tck50as,tck50rs = RunPrep(flagD,Vn)

    if len(species_initializations)==0:
        species_sheet = np.array([np.array(line.strip().split("\t")) for line in open('Species.txt', encoding='latin-1')])
        species_initializations = []
        for row in species_sheet[1:]:
            species_initializations.append(float(row[2]))
        species_initializations = np.array(species_initializations)

    xoutS_all = np.zeros(shape=(NSteps+1,len(species_initializations)))
    xoutS_all[0,:] = species_initializations # 24hr time point

    if len(genedata)==0:
        genedata = genedata0
    xoutG_all = np.zeros(shape=(NSteps+1,len(genedata)))
    xoutG_all[0,:] = genedata

    solver = model.getSolver() # Create solver instance
    solver.setMaxSteps = 1e10

    for qq in range(NSteps):
        genedata,xmNnM,AllGenesVec = SGEmodule(flagD,ts,xoutG_all[qq,:],xoutS_all[qq,:],Vn,Vc,kTCmaxs,kTCleak,kTCd,AllGenesVec, \
                                        GenePositionMatrix,TAs0,TRs0,kGin_1,kGac_1,tcnas,tck50as,tcnrs,tck50rs)
        xoutS_all[qq,773:914] = xmNnM # Update mRNA concentrations
        model.setInitialStates(xoutS_all[qq,:]) # Set new starting ICs
        rdata = amici.runAmiciSimulation(model, solver)  # Run simulation
        xoutS_all[qq+1,:] = rdata['x'][-1,:]
        xoutG_all[qq+1,:] = genedata
        if rdata['x'][-1,103] < rdata['x'][-1,105]:
            print('Apoptosis happened')
            break
    xoutS_all = xoutS_all[~np.all(xoutS_all == 0, axis=1)]
    xoutG_all = xoutG_all[~np.all(xoutG_all == 0, axis=1)]
    tout_all = tout_all[0:len(xoutS_all)]

    return xoutS_all, xoutG_all, tout_all
