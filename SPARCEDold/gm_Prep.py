import numpy as np
import scipy.stats
from random import *
import pandas as pd

def gm_Prep(flagD, gExp_mpc, mExp_mpc, kGin, kGac, kTCleak, kTCmaxs, kTCd):
    # Inputs:
    # flagD = Determistic(1) or Stochastic(0) simulation
    # gExp_mpc = molecules per cell gene numbers
    # mExp_mpc = molecules per cell mRNA numbers
    # kGin = rate of gene inactivation
    # kGac = rate of gene activation
    # kTCleak = Transcription leakage rates
    # kTCmaxs = Maximal transcription rates
    
    # Outputs:
    # genedata = active genes, inactive genes, mRNA numbers
    # GenePositionMatrix = matrix showing which genes are active
    # AllGenesVec = Number of active genes
    # xgac_mpc_D = molecules per cell active gene numbers - deterministic
    # xgin_mpc_D = molecules per cell inactive gene numbers - deterministic
    # xgac_mpc = molecules per cell active gene numbers - stochastic
    # xgin_mpc = molecules per cell inactive gene numbers - stochastic
    # kTCleak2 = Transcription leakage rates

    numberofgenes = int(len(gExp_mpc))
    a = np.arange(5,9)
    b = np.arange(12,25)
    indsD = np.concatenate((a, b), axis=None)
    mExp_mpc[indsD] = 17.0 # modify cell cycle gene mRNA numbers to 17

    ss = int(sum(gExp_mpc))
    GenePositionMatrix = np.zeros((numberofgenes,ss))
    ind = 0
    for i in range(numberofgenes):
        GenePositionMatrix[i, ind:ind+int(gExp_mpc[i])] = 1.0
        ind = ind + int(gExp_mpc[i])

    # xg Deterministic
    xgac_mpc_D = (kGac*gExp_mpc)/(kGin+kGac) #active genes initial condition
    xgin_mpc_D = gExp_mpc - xgac_mpc_D #inactive genes initial condition

    # xg Stochastic
    AllGenesVec = np.zeros(shape = (ss,1))
    IndsGenesOn = np.random.choice(ss, size=int(round(ss*kGac[0]/kGin[0])), replace=False)
    AllGenesVec[IndsGenesOn] = 1.0

    # Calculate Concentration of Active and Inactive Genes for each gene
    xgac_mpc = np.dot(GenePositionMatrix,AllGenesVec)
    xgac_mpc = xgac_mpc.ravel()
    xgin_mpc = gExp_mpc-xgac_mpc
    xgin_mpc = xgin_mpc.ravel()
    
    # kTCleak (deterministic)
    aa = np.multiply(kTCd,mExp_mpc)
    kTCleak2 = np.divide(aa,xgac_mpc_D)
    kTCleak2[np.isnan(kTCleak2)] = 0.0
    kTCleak2[np.isinf(kTCleak2)] = 0.0

    genedata = []
    if flagD==1:
        genedata = np.concatenate((xgac_mpc_D, xgin_mpc_D), axis=None)
        genedata = np.concatenate((genedata, mExp_mpc), axis=None)
    else:
        genedata = np.concatenate((xgac_mpc, xgin_mpc), axis=None)
        genedata = np.concatenate((genedata, mExp_mpc), axis=None)
    
    return genedata, GenePositionMatrix, AllGenesVec, xgac_mpc_D, xgin_mpc_D, xgac_mpc, xgin_mpc, kTCleak2