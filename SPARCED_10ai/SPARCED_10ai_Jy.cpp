#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "SPARCED_10ai_k.h"
#include "SPARCED_10ai_y.h"
#include "SPARCED_10ai_sigmay.h"
#include "SPARCED_10ai_my.h"

namespace amici {
namespace model_SPARCED_10ai {

void Jy_SPARCED_10ai(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_p53, 2)) + 0.5*std::pow(-mp53 + p53, 2)/std::pow(sigma_p53, 2);
            break;
        case 1:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_Mdm2, 2)) + 0.5*std::pow(Mdm2 - mMdm2, 2)/std::pow(sigma_Mdm2, 2);
            break;
        case 2:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_Wip1, 2)) + 0.5*std::pow(Wip1 - mWip1, 2)/std::pow(sigma_Wip1, 2);
            break;
        case 3:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_BRCA2, 2)) + 0.5*std::pow(BRCA2 - mBRCA2, 2)/std::pow(sigma_BRCA2, 2);
            break;
        case 4:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_MSH6, 2)) + 0.5*std::pow(MSH6 - mMSH6, 2)/std::pow(sigma_MSH6, 2);
            break;
        case 5:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_MGMT, 2)) + 0.5*std::pow(MGMT - mMGMT, 2)/std::pow(sigma_MGMT, 2);
            break;
        case 6:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ARF, 2)) + 0.5*std::pow(ARF - mARF, 2)/std::pow(sigma_ARF, 2);
            break;
        case 7:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_MDM4, 2)) + 0.5*std::pow(MDM4 - mMDM4, 2)/std::pow(sigma_MDM4, 2);
            break;
        case 8:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ATM, 2)) + 0.5*std::pow(ATM - mATM, 2)/std::pow(sigma_ATM, 2);
            break;
        case 9:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ATR, 2)) + 0.5*std::pow(ATR - mATR, 2)/std::pow(sigma_ATR, 2);
            break;
        case 10:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_RB, 2)) + 0.5*std::pow(RB - mRB, 2)/std::pow(sigma_RB, 2);
            break;
        case 11:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_E2F, 2)) + 0.5*std::pow(E2F - mE2F, 2)/std::pow(sigma_E2F, 2);
            break;
        case 12:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_Cd, 2)) + 0.5*std::pow(Cd - mCd, 2)/std::pow(sigma_Cd, 2);
            break;
        case 13:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_Ce, 2)) + 0.5*std::pow(Ce - mCe, 2)/std::pow(sigma_Ce, 2);
            break;
        case 14:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_Skp2, 2)) + 0.5*std::pow(Skp2 - mSkp2, 2)/std::pow(sigma_Skp2, 2);
            break;
        case 15:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_Pai, 2)) + 0.5*std::pow(Pai - mPai, 2)/std::pow(sigma_Pai, 2);
            break;
        case 16:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_Pei, 2)) + 0.5*std::pow(Pei - mPei, 2)/std::pow(sigma_Pei, 2);
            break;
        case 17:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_Pbi, 2)) + 0.5*std::pow(Pbi - mPbi, 2)/std::pow(sigma_Pbi, 2);
            break;
        case 18:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_Ca, 2)) + 0.5*std::pow(Ca - mCa, 2)/std::pow(sigma_Ca, 2);
            break;
        case 19:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_p27, 2)) + 0.5*std::pow(-mp27 + p27, 2)/std::pow(sigma_p27, 2);
            break;
        case 20:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_Cdh1a, 2)) + 0.5*std::pow(Cdh1a - mCdh1a, 2)/std::pow(sigma_Cdh1a, 2);
            break;
        case 21:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_Cb, 2)) + 0.5*std::pow(Cb - mCb, 2)/std::pow(sigma_Cb, 2);
            break;
        case 22:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_Cdc20, 2)) + 0.5*std::pow(Cdc20 - mCdc20, 2)/std::pow(sigma_Cdc20, 2);
            break;
        case 23:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_Wee1, 2)) + 0.5*std::pow(Wee1 - mWee1, 2)/std::pow(sigma_Wee1, 2);
            break;
        case 24:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_Chk1, 2)) + 0.5*std::pow(Chk1 - mChk1, 2)/std::pow(sigma_Chk1, 2);
            break;
        case 25:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_p21, 2)) + 0.5*std::pow(-mp21 + p21, 2)/std::pow(sigma_p21, 2);
            break;
        case 26:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_L, 2)) + 0.5*std::pow(L - mL, 2)/std::pow(sigma_L, 2);
            break;
        case 27:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_R, 2)) + 0.5*std::pow(R - mR, 2)/std::pow(sigma_R, 2);
            break;
        case 28:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_flip, 2)) + 0.5*std::pow(flip - mflip, 2)/std::pow(sigma_flip, 2);
            break;
        case 29:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_C8, 2)) + 0.5*std::pow(C8 - mC8, 2)/std::pow(sigma_C8, 2);
            break;
        case 30:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_Bar, 2)) + 0.5*std::pow(Bar - mBar, 2)/std::pow(sigma_Bar, 2);
            break;
        case 31:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_C3, 2)) + 0.5*std::pow(C3 - mC3, 2)/std::pow(sigma_C3, 2);
            break;
        case 32:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_C6, 2)) + 0.5*std::pow(C6 - mC6, 2)/std::pow(sigma_C6, 2);
            break;
        case 33:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_XIAP, 2)) + 0.5*std::pow(XIAP - mXIAP, 2)/std::pow(sigma_XIAP, 2);
            break;
        case 34:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_PARP, 2)) + 0.5*std::pow(PARP - mPARP, 2)/std::pow(sigma_PARP, 2);
            break;
        case 35:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_Bid, 2)) + 0.5*std::pow(Bid - mBid, 2)/std::pow(sigma_Bid, 2);
            break;
        case 36:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_Bcl2, 2)) + 0.5*std::pow(Bcl2 - mBcl2, 2)/std::pow(sigma_Bcl2, 2);
            break;
        case 37:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_Bax, 2)) + 0.5*std::pow(Bax - mBax, 2)/std::pow(sigma_Bax, 2);
            break;
        case 38:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_CytoC, 2)) + 0.5*std::pow(CytoC - mCytoC, 2)/std::pow(sigma_CytoC, 2);
            break;
        case 39:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_Smac, 2)) + 0.5*std::pow(Smac - mSmac, 2)/std::pow(sigma_Smac, 2);
            break;
        case 40:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_Apaf, 2)) + 0.5*std::pow(Apaf - mApaf, 2)/std::pow(sigma_Apaf, 2);
            break;
        case 41:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_C9, 2)) + 0.5*std::pow(C9 - mC9, 2)/std::pow(sigma_C9, 2);
            break;
        case 42:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_BAD, 2)) + 0.5*std::pow(BAD - mBAD, 2)/std::pow(sigma_BAD, 2);
            break;
        case 43:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_PUMA, 2)) + 0.5*std::pow(PUMA - mPUMA, 2)/std::pow(sigma_PUMA, 2);
            break;
        case 44:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_NOXA, 2)) + 0.5*std::pow(NOXA - mNOXA, 2)/std::pow(sigma_NOXA, 2);
            break;
        case 45:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_BIM, 2)) + 0.5*std::pow(BIM - mBIM, 2)/std::pow(sigma_BIM, 2);
            break;
        case 46:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_E, 2)) + 0.5*std::pow(E - mE, 2)/std::pow(sigma_E, 2);
            break;
        case 47:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_H, 2)) + 0.5*std::pow(H - mH, 2)/std::pow(sigma_H, 2);
            break;
        case 48:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_HGF, 2)) + 0.5*std::pow(HGF - mHGF, 2)/std::pow(sigma_HGF, 2);
            break;
        case 49:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_P, 2)) + 0.5*std::pow(P - mP, 2)/std::pow(sigma_P, 2);
            break;
        case 50:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_F, 2)) + 0.5*std::pow(F - mF, 2)/std::pow(sigma_F, 2);
            break;
        case 51:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_I, 2)) + 0.5*std::pow(I - mI, 2)/std::pow(sigma_I, 2);
            break;
        case 52:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_IN, 2)) + 0.5*std::pow(IN - mIN, 2)/std::pow(sigma_IN, 2);
            break;
        case 53:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_E1, 2)) + 0.5*std::pow(E1 - mE1, 2)/std::pow(sigma_E1, 2);
            break;
        case 54:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_E2, 2)) + 0.5*std::pow(E2 - mE2, 2)/std::pow(sigma_E2, 2);
            break;
        case 55:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_E3, 2)) + 0.5*std::pow(E3 - mE3, 2)/std::pow(sigma_E3, 2);
            break;
        case 56:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_E4, 2)) + 0.5*std::pow(E4 - mE4, 2)/std::pow(sigma_E4, 2);
            break;
        case 57:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_Ev3, 2)) + 0.5*std::pow(Ev3 - mEv3, 2)/std::pow(sigma_Ev3, 2);
            break;
        case 58:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_Met, 2)) + 0.5*std::pow(Met - mMet, 2)/std::pow(sigma_Met, 2);
            break;
        case 59:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_Pr, 2)) + 0.5*std::pow(Pr - mPr, 2)/std::pow(sigma_Pr, 2);
            break;
        case 60:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_Fr, 2)) + 0.5*std::pow(Fr - mFr, 2)/std::pow(sigma_Fr, 2);
            break;
        case 61:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_Ir, 2)) + 0.5*std::pow(Ir - mIr, 2)/std::pow(sigma_Ir, 2);
            break;
        case 62:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_Isr, 2)) + 0.5*std::pow(Isr - mIsr, 2)/std::pow(sigma_Isr, 2);
            break;
        case 63:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_IRS, 2)) + 0.5*std::pow(IRS - mIRS, 2)/std::pow(sigma_IRS, 2);
            break;
        case 64:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_Sp, 2)) + 0.5*std::pow(Sp - mSp, 2)/std::pow(sigma_Sp, 2);
            break;
        case 65:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_Cbl, 2)) + 0.5*std::pow(Cbl - mCbl, 2)/std::pow(sigma_Cbl, 2);
            break;
        case 66:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_G2, 2)) + 0.5*std::pow(G2 - mG2, 2)/std::pow(sigma_G2, 2);
            break;
        case 67:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_PLCg, 2)) + 0.5*std::pow(PLCg - mPLCg, 2)/std::pow(sigma_PLCg, 2);
            break;
        case 68:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_PI3KC1, 2)) + 0.5*std::pow(PI3KC1 - mPI3KC1, 2)/std::pow(sigma_PI3KC1, 2);
            break;
        case 69:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_PI3KR1, 2)) + 0.5*std::pow(PI3KR1 - mPI3KR1, 2)/std::pow(sigma_PI3KR1, 2);
            break;
        case 70:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_PI3K2, 2)) + 0.5*std::pow(PI3K2 - mPI3K2, 2)/std::pow(sigma_PI3K2, 2);
            break;
        case 71:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_GRP, 2)) + 0.5*std::pow(GRP - mGRP, 2)/std::pow(sigma_GRP, 2);
            break;
        case 72:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_Ras, 2)) + 0.5*std::pow(Ras - mRas, 2)/std::pow(sigma_Ras, 2);
            break;
        case 73:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_NF1, 2)) + 0.5*std::pow(NF1 - mNF1, 2)/std::pow(sigma_NF1, 2);
            break;
        case 74:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_CRaf, 2)) + 0.5*std::pow(CRaf - mCRaf, 2)/std::pow(sigma_CRaf, 2);
            break;
        case 75:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_BRaf, 2)) + 0.5*std::pow(BRaf - mBRaf, 2)/std::pow(sigma_BRaf, 2);
            break;
        case 76:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_MEK, 2)) + 0.5*std::pow(MEK - mMEK, 2)/std::pow(sigma_MEK, 2);
            break;
        case 77:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_MKP3, 2)) + 0.5*std::pow(MKP3 - mMKP3, 2)/std::pow(sigma_MKP3, 2);
            break;
        case 78:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_RSK, 2)) + 0.5*std::pow(RSK - mRSK, 2)/std::pow(sigma_RSK, 2);
            break;
        case 79:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_MKP1, 2)) + 0.5*std::pow(MKP1 - mMKP1, 2)/std::pow(sigma_MKP1, 2);
            break;
        case 80:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_Fos, 2)) + 0.5*std::pow(Fos - mFos, 2)/std::pow(sigma_Fos, 2);
            break;
        case 81:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_Jun, 2)) + 0.5*std::pow(Jun - mJun, 2)/std::pow(sigma_Jun, 2);
            break;
        case 82:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_Myc, 2)) + 0.5*std::pow(Myc - mMyc, 2)/std::pow(sigma_Myc, 2);
            break;
        case 83:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_bCATENIN, 2)) + 0.5*std::pow(bCATENIN - mbCATENIN, 2)/std::pow(sigma_bCATENIN, 2);
            break;
        case 84:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_PTEN, 2)) + 0.5*std::pow(PTEN - mPTEN, 2)/std::pow(sigma_PTEN, 2);
            break;
        case 85:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_AKT, 2)) + 0.5*std::pow(AKT - mAKT, 2)/std::pow(sigma_AKT, 2);
            break;
        case 86:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_PDK1, 2)) + 0.5*std::pow(PDK1 - mPDK1, 2)/std::pow(sigma_PDK1, 2);
            break;
        case 87:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_Rictor, 2)) + 0.5*std::pow(Rictor - mRictor, 2)/std::pow(sigma_Rictor, 2);
            break;
        case 88:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_mTOR, 2)) + 0.5*std::pow(mTOR - mmTOR, 2)/std::pow(sigma_mTOR, 2);
            break;
        case 89:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_GSK3b, 2)) + 0.5*std::pow(GSK3b - mGSK3b, 2)/std::pow(sigma_GSK3b, 2);
            break;
        case 90:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_TSC1, 2)) + 0.5*std::pow(TSC1 - mTSC1, 2)/std::pow(sigma_TSC1, 2);
            break;
        case 91:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_TSC2, 2)) + 0.5*std::pow(TSC2 - mTSC2, 2)/std::pow(sigma_TSC2, 2);
            break;
        case 92:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_PKC, 2)) + 0.5*std::pow(PKC - mPKC, 2)/std::pow(sigma_PKC, 2);
            break;
        case 93:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_RKIP, 2)) + 0.5*std::pow(RKIP - mRKIP, 2)/std::pow(sigma_RKIP, 2);
            break;
        case 94:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERK, 2)) + 0.5*std::pow(ERK - mERK, 2)/std::pow(sigma_ERK, 2);
            break;
        case 95:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_FOXO, 2)) + 0.5*std::pow(FOXO - mFOXO, 2)/std::pow(sigma_FOXO, 2);
            break;
        case 96:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_Rheb, 2)) + 0.5*std::pow(Rheb - mRheb, 2)/std::pow(sigma_Rheb, 2);
            break;
        case 97:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_Raptor, 2)) + 0.5*std::pow(Raptor - mRaptor, 2)/std::pow(sigma_Raptor, 2);
            break;
        case 98:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S6K, 2)) + 0.5*std::pow(S6K - mS6K, 2)/std::pow(sigma_S6K, 2);
            break;
        case 99:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_EIF4EBP1, 2)) + 0.5*std::pow(EIF4EBP1 - mEIF4EBP1, 2)/std::pow(sigma_EIF4EBP1, 2);
            break;
        case 100:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_SOS, 2)) + 0.5*std::pow(SOS - mSOS, 2)/std::pow(sigma_SOS, 2);
            break;
        case 101:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_EIF4E, 2)) + 0.5*std::pow(EIF4E - mEIF4E, 2)/std::pow(sigma_EIF4E, 2);
            break;
    }
}

} // namespace model_SPARCED_10ai
} // namespace amici
