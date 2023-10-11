import numpy as np
import scipy.stats
from random import *
import pandas as pd

def SGEmodule(flagD,ts,genedata,spdata,Vn,Vc,kTCmaxs,kTCleak,kTCd,AllGenesVec,GenePositionMatrix,kGin_1,kGac_1, 
              tcnas,tck50as,tcnrs,tck50rs,spIDs,mRNAInds0):
    # Inputs:
    # flagD = deterministic (1) or stochastic (0) simulation
    # ts = time
    # genedata = 1D array of latest gene-expression module concentrations
    # spdata = 1D array of latest species concentrations
    # Vn = nuclear volume
    # Vc = cytoplasmic volume

    # Outputs:
    # genedataNew = new active genes, inactive genes, and mRNA concentrations
    # AllGenesVecNew = New array of all genes
    # kmRNAlist = list of names for paramters to change with new mRNA concentrations

    mpc2nmcf_Vn = 1.0E9/(Vn*6.023E+23)
    mpc2nmcf_Vc = 1.0E9/(Vc*6.023E+23)
    numberofgenes = len(tcnas)
    numberofTARs = len(tcnas[0])
    
    ix = 0
    xgac = genedata[ix:ix+numberofgenes] # active genes
    ix = ix+numberofgenes
    xgin = genedata[ix:ix+numberofgenes] # inactive genes
    xm = np.divide(spdata[mRNAInds0:],mpc2nmcf_Vc) # mRNAs: nM --> molecules per cell
        
    TARarr = np.array(spdata[spIDs])
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
    
    # vTC
    hills = np.matrix(hills)
    induced = np.multiply(np.multiply(xgac,kTCmaxs),hills)
    induced = induced.flatten()
    leak = np.multiply(xgac,kTCleak)
    vTC = np.add(leak,induced)
    vTC.flatten()
    vTC = np.squeeze(np.asarray(vTC))

    # vTCd
    vTCd= np.transpose(np.multiply(kTCd,xm));
    vTCd = np.squeeze(np.asarray(vTCd))
    
    # If deterministic simulation:
    if flagD: 
        Nb = vTC*ts;
        Nd = vTCd*ts;
        xgacN = genedata[0:numberofgenes]
        xginN = genedata[numberofgenes:numberofgenes*2]
        AllGenesVecN = []
    else:
        # Poisson processes
        poff = scipy.stats.poisson.pmf(0,kGin_1*ts)
        pon = scipy.stats.poisson.pmf(0,kGac_1*ts)

        # Generating random numbers and deciding which genes should turn off and on
        RandomNumbers = np.random.uniform(0,1,len(AllGenesVec))
        geneson = AllGenesVec.astype(bool).astype(int)
        genesoff = np.logical_not(geneson).astype(int)
        ac2in = np.logical_and(np.transpose(geneson.flatten()),RandomNumbers>=poff)
        in2ac = np.logical_and(np.transpose(genesoff.flatten()),RandomNumbers>=pon)

        # Allocate active and inactive genes
        AllGenesVecN = AllGenesVec
        AllGenesVecN[ac2in] = 0.0
        AllGenesVecN[in2ac] = 1.0

        xgacN = np.dot(GenePositionMatrix,AllGenesVecN)
        xgacN = xgacN.ravel()
        xginN = np.subtract(np.add(xgac,xgin),xgacN.flatten())
        xginN = xginN.ravel()

        # mRNA
        Nb = np.random.poisson(np.float64(vTC*ts))
        Nb = Nb.ravel()
        Nd = np.random.poisson(np.float64(vTCd*ts))
        Nd = Nd.ravel()
        # These genes and mRNAs we do not allow to fluctuate
        indsD = np.array([5,6,7,8,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29])
        Nb[indsD] = vTC[indsD]*ts
        Nd[indsD] = vTCd[indsD]*ts
        xgacN[indsD] = genedata[indsD]
        xginN[indsD] = genedata[indsD+numberofgenes]

    # Finish mRNA
    xmN = xm+Nb-Nd
    xmN[xmN<0.0] = 0.0
    xmN_nM = xmN*mpc2nmcf_Vc

    genedataNew = []
    genedataNew = np.concatenate((xgacN, xginN), axis=None)

    return genedataNew, xmN, AllGenesVecN
