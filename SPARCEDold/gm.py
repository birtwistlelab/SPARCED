import numpy as np
import scipy.stats
from random import *
import pandas as pd

def gm(flagD,ts,genedata,spdata,Vn,Vc,kGin,kGac,kTCmaxs,kTCleak,kTCd,AllGenesVec,GenePositionMatrix,TAs0,TRs0):
    # Inputs:
    # flagD = deterministic (1) or stochastic (0) simulation
    # ts = time step
    # genedata = 1D array of latest gene-expression module concentrations
    # spdata = 1D array of latest species concentrations
    # Vn = nuclear volume
    # Vc = cytoplasmic volume
    # kGin = rate of gene inactivation
    # kGac = rate of gene activation
    # kTCmaxs = Maximal transcription rates
    # kTCleak = Transcription leakage rates
    # kTCd = mRNA degradation rates
    # AllGenesVec = Number of active genes
    # GenePositionMatrix = matrix showing which genes are active
    # TAs0 = Genes x Transcriptional activators
    # TRs0 = Genes x Transcriptional repressors
   
    # Outputs:
    # genedataNew = new active genes, inactive genes, and mRNA concentrations
    # AllGenesVecNew = New array of all genes
    
    mpc2nmcf_Vn = 1E9/(Vn*6.023E+23)
    numberofgenes = int(len(kTCleak))
    numberofTARs = 7
    a = np.arange(5,9)
    b = np.arange(12,30)
    indsD = np.concatenate((a, b), axis=None)

    # gm species
    a = 0
    xgac = genedata[a:a+numberofgenes]
    a = a+numberofgenes
    xgin = genedata[a:a+numberofgenes]
    a = a+numberofgenes
    xm = genedata[a:a+numberofgenes]

    # Gene switching constants
    kGin_1 = kGin[0]
    kGac_1 = kGac[0]

    # Get the latest concentrations of TARs
    pcFos_cJun = spdata[684] #1
    cMyc = spdata[685] #2
    p53ac = spdata[2] #3
    FOXOnuc = spdata[767] #4
    ppERKnuc = spdata[675] #5
    pRSKnuc = spdata[678] #6
    bCATENINnuc = spdata[686] #7

    cc = np.multiply(TAs0,np.array([pcFos_cJun, cMyc, p53ac, FOXOnuc, ppERKnuc, pRSKnuc, bCATENINnuc]))
    dd = cc*(1/mpc2nmcf_Vn) # convert to mpc from nM
    TAs = dd
    TAs.flatten()
    ee = np.multiply(TRs0,np.array([pcFos_cJun, 1, 1, 1, 1, 1, 1]))
    ff = ee*(1/mpc2nmcf_Vn)   
    TRs = ff
    TRs.flatten()
    
    tcnas = TAs0.copy() # used in Hill equations
    tcnas += 1.0
    tck50as = TAs0.copy()
    tcnrs = TRs0.copy()
    tcnrs += 1.0
    tck50rs = TRs0.copy()
    
    tcnas[[val for val in np.nonzero(tcnas[:,0]==2)],0] = 3.0 # pcFos_cJun
    tcnas[[val for val in np.nonzero(tcnas[:,1]==2)],1] = 3.0 # cMyc
    tcnas[[val for val in np.nonzero(tcnas[:,2]==2)],2] = 4.0 # p53ac
    tcnas[[val for val in np.nonzero(tcnas[:,3]==2)],3] = 4.0 # FOXOnuc
    tcnas[[val for val in np.nonzero(tcnas[:,4]==2)],4] = 4.0 # ppERKnuc
    tcnas[[val for val in np.nonzero(tcnas[:,5]==2)],5] = 4.0 # pRSKnuc
    tcnas[[val for val in np.nonzero(tcnas[:,6]==2)],6] = 4.0 # bCATENIN

    tcnrs[[val for val in np.nonzero(tcnrs[:,0]==2)],0] = 4.0 # pcFos_cJun

    # pcFos_cJun
    tck50as[9:12,0] = 1.25 # CyclinD
    tck50as[98,0] = 0.8 # cJun
    # cMyc
    tck50as[9:12,1] = 450.0 # CyclinD
    # p53ac
    tck50as[25,2] = 50.0 # p21
    tck50as[52:54,2] = 1350.0 # PUMA,NOXA
    # FOXO
    tck50as[54,3] = 45.0 # 19
    tck50as[[57,58,59,60,62,64,65,126,127,135,139],3] = 60.0 # RTKs
    # ppERKnuc
    tck50as[67,4] = 65.0 # SPRY2
    tck50as[[91,96],4] = 40.0 # DUSPs
    tck50as[97,4] = 20.0 # cFos
    # pRSKnuc
    tck50as[67,5] = 20.0 # SPRY2
    tck50as[[91,96],5] = 10.0 # DUSPs
    tck50as[97,5] = 5.0 # cFos
    # bCATENINnuc
    tck50as[99,6] = 250.0 # cMyc

    tck50rs[97,0] = tck50as[98,0] # cFos

    # Convert to molecules per cell
    tck50as = tck50as*(1/mpc2nmcf_Vn)
    tck50rs = tck50rs*(1/mpc2nmcf_Vn)

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
    # hills = np.matrix.transpose(hills)
    induced = np.multiply(np.multiply(xgac,kTCmaxs),hills)
    induced = induced.flatten()
    leak = np.multiply(xgac,kTCleak)
    vTC = np.add(leak,induced)
    vTC.flatten()
    vTC = np.squeeze(np.asarray(vTC))

    # vTCd
    vTCd= np.transpose(np.multiply(kTCd,xm));
    vTCd = np.squeeze(np.asarray(vTCd))

    # Poisson Stuff
    poff = scipy.stats.poisson.pmf(0,kGin_1*ts)
    pon = scipy.stats.poisson.pmf(0,kGac_1*ts)

    # Generating random numbers and deciding which genes should turn off and on
    RandomNumbers = np.random.uniform(0,1,len(AllGenesVec))
    geneson = AllGenesVec.astype(bool).astype(int)
    genesoff = np.logical_not(geneson).astype(int)
    ac2in = np.logical_and(np.transpose(geneson.flatten()),RandomNumbers>=poff)
    in2ac = np.logical_and(np.transpose(genesoff.flatten()),RandomNumbers>=pon)

    # Generating new AllGenesVec and allocating active and inactive genes
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
    # These genes and mRNAs we don't allow to fluctuate
    Nb[indsD] = vTC[indsD]*ts
    Nd[indsD] = vTCd[indsD]*ts
    xgacN[indsD] = genedata[indsD]
    xginN[indsD] = genedata[indsD+numberofgenes]

    # If deterministic simulation:
    if flagD: 
        Nb = vTC*ts;
        Nd = vTCd*ts;
        xgacN = genedata[0:numberofgenes]
        xginN = genedata[numberofgenes:numberofgenes*2]

    # Finish mRNA
    xmN = xm+Nb-Nd
    xmN[xmN<0.0] = 0.0
    xmN_nM = xmN*(1E9/(Vc*6.023E+23))

    genedataNew = []
    genedataNew = np.concatenate((xgacN, xginN), axis=None)
    genedataNew = np.concatenate((genedataNew, xmN), axis=None)

    return genedataNew, AllGenesVecN