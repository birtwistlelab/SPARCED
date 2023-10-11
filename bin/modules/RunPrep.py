import numpy as np
import pandas as pd
import os

def RunPrep(flagD,Vn,model):
    
    cd = os.getcwd()
    wd = os.path.dirname(cd)
    input_path = os.path.join(wd,'input_files')
    
    
    kGsRead = pd.read_csv(os.path.join(input_path,'OmicsData.txt'),header=0,index_col=0,sep="\t")
    gExp_mpc = np.float64(kGsRead.values[:,0])
    mExp_mpc = np.float64(kGsRead.values[:,1])
    kGin = np.float64(kGsRead.values[:,2])
    kGac = np.float64(kGsRead.values[:,3])
    kTCleak = np.float64(kGsRead.values[:,4])
    kTCmaxs = np.float64(kGsRead.values[:,5])
    kTCd = np.float64(kGsRead.values[:,6])

    # Read-in the activators matrix and assign concentrations of activators
    TARsRead = pd.read_csv(os.path.join(input_path,'GeneReg.txt'),header=0,index_col=0,sep="\t")
    TARs0 = (TARsRead.values)
    numberofTARs = len(TARsRead.columns)
    spnames = [ele for ele in model.getStateIds()]
    spIDs = []
    for qq in range(numberofTARs):
        sps = spnames.index(TARsRead.columns[qq]) 
        spIDs.append(sps)
    TARsRead = None
    
    numberofgenes = int(len(gExp_mpc))
    indsDm = [5,6,7,8,12,13,14,15,16,17,18,19,20,21,22,23,24]
    mExp_mpc[indsDm] = 17.0 # modify cell cycle gene mRNA numbers to 17

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

    genedata = []
    if flagD==1:
        genedata = np.concatenate((xgac_mpc_D, xgin_mpc_D), axis=None)
    else:
        genedata = np.concatenate((xgac_mpc, xgin_mpc), axis=None)
    
    # Gene switching constants
    kGin_1 = kGin[0]
    kGac_1 = kGac[0]
    
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
    
    return genedata, GenePositionMatrix, AllGenesVec, kTCmaxs, kTCleak, kGin_1, kGac_1, kTCd, TARs0, tcnas, tcnrs, tck50as, tck50rs, spIDs 
