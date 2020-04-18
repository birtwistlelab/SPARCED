import numpy as np
import pandas as pd

def RunPrep(flagD,Vn, input_data_folder):
    # Inputs:
    # flagD = deterministic (1) or stochastic (0) simulation
    # Vn = nuclear volume

    # Outputs:
    # genedata = 1D array of latest gene-expression module concentrations
    # mExp_mpc = molecules per cell mRNA numbers
    # GenePositionMatrix = matrix showing which genes are active
    # AllGenesVec = Number of active genes
    # kTCmaxs = Maximal transcription rates
    # kTCleak = Transcription leakage rates - from data
    # kTCleak2 = Transcription leakage rates - calculated
    # kGin_1 = rate of gene inactivation
    # kGac_1 = rate of gene activation
    # kTCd = mRNA degradation rates
    # TAs_data = Genes x Transcriptional activators
    # TRs_data = Genes x Transcriptional repressors
    # tcnas = Hill coefficients for transcriptional activators
    # tck50as = K50 values for transcriptional activators
    # tcnrs = Hill coefficients for transcriptional repressors
    # tck50rs = K50 values for transcriptional repressors

    kGsRead_sheet = np.array([np.array(line.strip().split("\t")) for line in open(input_data_folder+'kGeneMod.txt', encoding='latin-1')])
    kGsRead_data = []
    for row in kGsRead_sheet[1:]:
        kGsRead_data.append(row[1:])
    kGsRead_data = np.array(kGsRead_data)

    gExp_mpc = np.float64(kGsRead_data[:,0])
    mExp_mpc = np.float64(kGsRead_data[:,1])
    kGin = np.float64(kGsRead_data[:,2])
    kGac = np.float64(kGsRead_data[:,3])
    kTCleak = np.float64(kGsRead_data[:,4])
    kTCmaxs = np.float64(kGsRead_data[:,5])
    kTCd = np.float64(kGsRead_data[:,6])

    # Read-in the activators matrix and assign concentrations of activators
    TAs_sheet = np.array([np.array(line.strip().split(",")) for line in open(input_data_folder+'TAs.csv', encoding='latin-1')])
    TAs_data = []
    for row in TAs_sheet[1:]:
        TAs_data.append(row[1:])
    TAs_data = np.array(TAs_data, dtype=float)


    # Read-in the repressors matrix and assign concentrations of repressors
    TRs_sheet = np.array([np.array(line.strip().split(",")) for line in open(input_data_folder+'TRs.csv', encoding='latin-1')])
    TRs_data = []
    for row in TRs_sheet[1:]:
        TRs_data.append(row[1:])
    TRs_data = np.array(TRs_data, dtype=float)


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

    # kTCleak (deterministic)
    aa = np.multiply(kTCd,mExp_mpc)
    kTCleak2 = np.divide(aa,xgac_mpc_D)
    kTCleak2[np.isnan(kTCleak2)] = 0.0
    kTCleak2[np.isinf(kTCleak2)] = 0.0

    genedata = []
    if flagD==1:
        genedata = np.concatenate((xgac_mpc_D, xgin_mpc_D), axis=None)
    else:
        genedata = np.concatenate((xgac_mpc, xgin_mpc), axis=None)

    # Gene switching constants
    kGin_1 = kGin[0]
    kGac_1 = kGac[0]

    tcnas = TAs_data.copy() # used in Hill equations
    tcnas += 1.0
    tck50as = TAs_data.copy()
    tcnrs = TRs_data.copy()
    tcnrs += 1.0
    tck50rs = TRs_data.copy()

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

    mpc2nmcf_Vn = 1.0E9/(Vn*6.023E+23)
    # Convert to molecules per cell
    tck50as = tck50as*(1/mpc2nmcf_Vn)
    tck50rs = tck50rs*(1/mpc2nmcf_Vn)

    return genedata, mExp_mpc, GenePositionMatrix, AllGenesVec, kTCmaxs, kTCleak, kTCleak2, kGin_1, kGac_1, kTCd, TAs_data, TRs_data, tcnas, tcnrs, tck50as, tck50rs
