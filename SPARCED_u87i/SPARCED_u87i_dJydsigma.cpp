#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "SPARCED_u87i_k.h"
#include "SPARCED_u87i_y.h"
#include "SPARCED_u87i_sigmay.h"
#include "SPARCED_u87i_my.h"

namespace amici {
namespace model_SPARCED_u87i {

void dJydsigma_SPARCED_u87i(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = 1.0/sigma_p53 - 1.0*std::pow(-mp53 + p53, 2)/std::pow(sigma_p53, 3);
            break;
        case 1:
            dJydsigma[1] = 1.0/sigma_Mdm2 - 1.0*std::pow(Mdm2 - mMdm2, 2)/std::pow(sigma_Mdm2, 3);
            break;
        case 2:
            dJydsigma[2] = 1.0/sigma_Wip1 - 1.0*std::pow(Wip1 - mWip1, 2)/std::pow(sigma_Wip1, 3);
            break;
        case 3:
            dJydsigma[3] = 1.0/sigma_BRCA2 - 1.0*std::pow(BRCA2 - mBRCA2, 2)/std::pow(sigma_BRCA2, 3);
            break;
        case 4:
            dJydsigma[4] = 1.0/sigma_MSH6 - 1.0*std::pow(MSH6 - mMSH6, 2)/std::pow(sigma_MSH6, 3);
            break;
        case 5:
            dJydsigma[5] = 1.0/sigma_MGMT - 1.0*std::pow(MGMT - mMGMT, 2)/std::pow(sigma_MGMT, 3);
            break;
        case 6:
            dJydsigma[6] = 1.0/sigma_ARF - 1.0*std::pow(ARF - mARF, 2)/std::pow(sigma_ARF, 3);
            break;
        case 7:
            dJydsigma[7] = 1.0/sigma_MDM4 - 1.0*std::pow(MDM4 - mMDM4, 2)/std::pow(sigma_MDM4, 3);
            break;
        case 8:
            dJydsigma[8] = 1.0/sigma_ATM - 1.0*std::pow(ATM - mATM, 2)/std::pow(sigma_ATM, 3);
            break;
        case 9:
            dJydsigma[9] = 1.0/sigma_ATR - 1.0*std::pow(ATR - mATR, 2)/std::pow(sigma_ATR, 3);
            break;
        case 10:
            dJydsigma[10] = 1.0/sigma_RB - 1.0*std::pow(RB - mRB, 2)/std::pow(sigma_RB, 3);
            break;
        case 11:
            dJydsigma[11] = 1.0/sigma_E2F - 1.0*std::pow(E2F - mE2F, 2)/std::pow(sigma_E2F, 3);
            break;
        case 12:
            dJydsigma[12] = 1.0/sigma_Cd - 1.0*std::pow(Cd - mCd, 2)/std::pow(sigma_Cd, 3);
            break;
        case 13:
            dJydsigma[13] = 1.0/sigma_Ce - 1.0*std::pow(Ce - mCe, 2)/std::pow(sigma_Ce, 3);
            break;
        case 14:
            dJydsigma[14] = 1.0/sigma_Skp2 - 1.0*std::pow(Skp2 - mSkp2, 2)/std::pow(sigma_Skp2, 3);
            break;
        case 15:
            dJydsigma[15] = 1.0/sigma_Pai - 1.0*std::pow(Pai - mPai, 2)/std::pow(sigma_Pai, 3);
            break;
        case 16:
            dJydsigma[16] = 1.0/sigma_Pei - 1.0*std::pow(Pei - mPei, 2)/std::pow(sigma_Pei, 3);
            break;
        case 17:
            dJydsigma[17] = 1.0/sigma_Pbi - 1.0*std::pow(Pbi - mPbi, 2)/std::pow(sigma_Pbi, 3);
            break;
        case 18:
            dJydsigma[18] = 1.0/sigma_Ca - 1.0*std::pow(Ca - mCa, 2)/std::pow(sigma_Ca, 3);
            break;
        case 19:
            dJydsigma[19] = 1.0/sigma_p27 - 1.0*std::pow(-mp27 + p27, 2)/std::pow(sigma_p27, 3);
            break;
        case 20:
            dJydsigma[20] = 1.0/sigma_Cdh1a - 1.0*std::pow(Cdh1a - mCdh1a, 2)/std::pow(sigma_Cdh1a, 3);
            break;
        case 21:
            dJydsigma[21] = 1.0/sigma_Cb - 1.0*std::pow(Cb - mCb, 2)/std::pow(sigma_Cb, 3);
            break;
        case 22:
            dJydsigma[22] = 1.0/sigma_Cdc20 - 1.0*std::pow(Cdc20 - mCdc20, 2)/std::pow(sigma_Cdc20, 3);
            break;
        case 23:
            dJydsigma[23] = 1.0/sigma_Wee1 - 1.0*std::pow(Wee1 - mWee1, 2)/std::pow(sigma_Wee1, 3);
            break;
        case 24:
            dJydsigma[24] = 1.0/sigma_Chk1 - 1.0*std::pow(Chk1 - mChk1, 2)/std::pow(sigma_Chk1, 3);
            break;
        case 25:
            dJydsigma[25] = 1.0/sigma_p21 - 1.0*std::pow(-mp21 + p21, 2)/std::pow(sigma_p21, 3);
            break;
        case 26:
            dJydsigma[26] = 1.0/sigma_L - 1.0*std::pow(L - mL, 2)/std::pow(sigma_L, 3);
            break;
        case 27:
            dJydsigma[27] = 1.0/sigma_R - 1.0*std::pow(R - mR, 2)/std::pow(sigma_R, 3);
            break;
        case 28:
            dJydsigma[28] = 1.0/sigma_flip - 1.0*std::pow(flip - mflip, 2)/std::pow(sigma_flip, 3);
            break;
        case 29:
            dJydsigma[29] = 1.0/sigma_C8 - 1.0*std::pow(C8 - mC8, 2)/std::pow(sigma_C8, 3);
            break;
        case 30:
            dJydsigma[30] = 1.0/sigma_Bar - 1.0*std::pow(Bar - mBar, 2)/std::pow(sigma_Bar, 3);
            break;
        case 31:
            dJydsigma[31] = 1.0/sigma_C3 - 1.0*std::pow(C3 - mC3, 2)/std::pow(sigma_C3, 3);
            break;
        case 32:
            dJydsigma[32] = 1.0/sigma_C6 - 1.0*std::pow(C6 - mC6, 2)/std::pow(sigma_C6, 3);
            break;
        case 33:
            dJydsigma[33] = 1.0/sigma_XIAP - 1.0*std::pow(XIAP - mXIAP, 2)/std::pow(sigma_XIAP, 3);
            break;
        case 34:
            dJydsigma[34] = 1.0/sigma_PARP - 1.0*std::pow(PARP - mPARP, 2)/std::pow(sigma_PARP, 3);
            break;
        case 35:
            dJydsigma[35] = 1.0/sigma_Bid - 1.0*std::pow(Bid - mBid, 2)/std::pow(sigma_Bid, 3);
            break;
        case 36:
            dJydsigma[36] = 1.0/sigma_Bcl2 - 1.0*std::pow(Bcl2 - mBcl2, 2)/std::pow(sigma_Bcl2, 3);
            break;
        case 37:
            dJydsigma[37] = 1.0/sigma_Bax - 1.0*std::pow(Bax - mBax, 2)/std::pow(sigma_Bax, 3);
            break;
        case 38:
            dJydsigma[38] = 1.0/sigma_CytoC - 1.0*std::pow(CytoC - mCytoC, 2)/std::pow(sigma_CytoC, 3);
            break;
        case 39:
            dJydsigma[39] = 1.0/sigma_Smac - 1.0*std::pow(Smac - mSmac, 2)/std::pow(sigma_Smac, 3);
            break;
        case 40:
            dJydsigma[40] = 1.0/sigma_Apaf - 1.0*std::pow(Apaf - mApaf, 2)/std::pow(sigma_Apaf, 3);
            break;
        case 41:
            dJydsigma[41] = 1.0/sigma_C9 - 1.0*std::pow(C9 - mC9, 2)/std::pow(sigma_C9, 3);
            break;
        case 42:
            dJydsigma[42] = 1.0/sigma_BAD - 1.0*std::pow(BAD - mBAD, 2)/std::pow(sigma_BAD, 3);
            break;
        case 43:
            dJydsigma[43] = 1.0/sigma_PUMA - 1.0*std::pow(PUMA - mPUMA, 2)/std::pow(sigma_PUMA, 3);
            break;
        case 44:
            dJydsigma[44] = 1.0/sigma_NOXA - 1.0*std::pow(NOXA - mNOXA, 2)/std::pow(sigma_NOXA, 3);
            break;
        case 45:
            dJydsigma[45] = 1.0/sigma_BIM - 1.0*std::pow(BIM - mBIM, 2)/std::pow(sigma_BIM, 3);
            break;
        case 46:
            dJydsigma[46] = 1.0/sigma_E - 1.0*std::pow(E - mE, 2)/std::pow(sigma_E, 3);
            break;
        case 47:
            dJydsigma[47] = 1.0/sigma_H - 1.0*std::pow(H - mH, 2)/std::pow(sigma_H, 3);
            break;
        case 48:
            dJydsigma[48] = 1.0/sigma_HGF - 1.0*std::pow(HGF - mHGF, 2)/std::pow(sigma_HGF, 3);
            break;
        case 49:
            dJydsigma[49] = 1.0/sigma_P - 1.0*std::pow(P - mP, 2)/std::pow(sigma_P, 3);
            break;
        case 50:
            dJydsigma[50] = 1.0/sigma_F - 1.0*std::pow(F - mF, 2)/std::pow(sigma_F, 3);
            break;
        case 51:
            dJydsigma[51] = 1.0/sigma_I - 1.0*std::pow(I - mI, 2)/std::pow(sigma_I, 3);
            break;
        case 52:
            dJydsigma[52] = 1.0/sigma_IN - 1.0*std::pow(IN - mIN, 2)/std::pow(sigma_IN, 3);
            break;
        case 53:
            dJydsigma[53] = 1.0/sigma_E1 - 1.0*std::pow(E1 - mE1, 2)/std::pow(sigma_E1, 3);
            break;
        case 54:
            dJydsigma[54] = 1.0/sigma_E2 - 1.0*std::pow(E2 - mE2, 2)/std::pow(sigma_E2, 3);
            break;
        case 55:
            dJydsigma[55] = 1.0/sigma_E3 - 1.0*std::pow(E3 - mE3, 2)/std::pow(sigma_E3, 3);
            break;
        case 56:
            dJydsigma[56] = 1.0/sigma_E4 - 1.0*std::pow(E4 - mE4, 2)/std::pow(sigma_E4, 3);
            break;
        case 57:
            dJydsigma[57] = 1.0/sigma_Ev3 - 1.0*std::pow(Ev3 - mEv3, 2)/std::pow(sigma_Ev3, 3);
            break;
        case 58:
            dJydsigma[58] = 1.0/sigma_Met - 1.0*std::pow(Met - mMet, 2)/std::pow(sigma_Met, 3);
            break;
        case 59:
            dJydsigma[59] = 1.0/sigma_Pr - 1.0*std::pow(Pr - mPr, 2)/std::pow(sigma_Pr, 3);
            break;
        case 60:
            dJydsigma[60] = 1.0/sigma_Fr - 1.0*std::pow(Fr - mFr, 2)/std::pow(sigma_Fr, 3);
            break;
        case 61:
            dJydsigma[61] = 1.0/sigma_Ir - 1.0*std::pow(Ir - mIr, 2)/std::pow(sigma_Ir, 3);
            break;
        case 62:
            dJydsigma[62] = 1.0/sigma_Isr - 1.0*std::pow(Isr - mIsr, 2)/std::pow(sigma_Isr, 3);
            break;
        case 63:
            dJydsigma[63] = 1.0/sigma_IRS - 1.0*std::pow(IRS - mIRS, 2)/std::pow(sigma_IRS, 3);
            break;
        case 64:
            dJydsigma[64] = 1.0/sigma_Sp - 1.0*std::pow(Sp - mSp, 2)/std::pow(sigma_Sp, 3);
            break;
        case 65:
            dJydsigma[65] = 1.0/sigma_Cbl - 1.0*std::pow(Cbl - mCbl, 2)/std::pow(sigma_Cbl, 3);
            break;
        case 66:
            dJydsigma[66] = 1.0/sigma_G2 - 1.0*std::pow(G2 - mG2, 2)/std::pow(sigma_G2, 3);
            break;
        case 67:
            dJydsigma[67] = 1.0/sigma_PLCg - 1.0*std::pow(PLCg - mPLCg, 2)/std::pow(sigma_PLCg, 3);
            break;
        case 68:
            dJydsigma[68] = 1.0/sigma_PI3KC1 - 1.0*std::pow(PI3KC1 - mPI3KC1, 2)/std::pow(sigma_PI3KC1, 3);
            break;
        case 69:
            dJydsigma[69] = 1.0/sigma_PI3KR1 - 1.0*std::pow(PI3KR1 - mPI3KR1, 2)/std::pow(sigma_PI3KR1, 3);
            break;
        case 70:
            dJydsigma[70] = 1.0/sigma_PI3K2 - 1.0*std::pow(PI3K2 - mPI3K2, 2)/std::pow(sigma_PI3K2, 3);
            break;
        case 71:
            dJydsigma[71] = 1.0/sigma_GRP - 1.0*std::pow(GRP - mGRP, 2)/std::pow(sigma_GRP, 3);
            break;
        case 72:
            dJydsigma[72] = 1.0/sigma_Ras - 1.0*std::pow(Ras - mRas, 2)/std::pow(sigma_Ras, 3);
            break;
        case 73:
            dJydsigma[73] = 1.0/sigma_NF1 - 1.0*std::pow(NF1 - mNF1, 2)/std::pow(sigma_NF1, 3);
            break;
        case 74:
            dJydsigma[74] = 1.0/sigma_CRaf - 1.0*std::pow(CRaf - mCRaf, 2)/std::pow(sigma_CRaf, 3);
            break;
        case 75:
            dJydsigma[75] = 1.0/sigma_BRaf - 1.0*std::pow(BRaf - mBRaf, 2)/std::pow(sigma_BRaf, 3);
            break;
        case 76:
            dJydsigma[76] = 1.0/sigma_MEK - 1.0*std::pow(MEK - mMEK, 2)/std::pow(sigma_MEK, 3);
            break;
        case 77:
            dJydsigma[77] = 1.0/sigma_MKP3 - 1.0*std::pow(MKP3 - mMKP3, 2)/std::pow(sigma_MKP3, 3);
            break;
        case 78:
            dJydsigma[78] = 1.0/sigma_RSK - 1.0*std::pow(RSK - mRSK, 2)/std::pow(sigma_RSK, 3);
            break;
        case 79:
            dJydsigma[79] = 1.0/sigma_MKP1 - 1.0*std::pow(MKP1 - mMKP1, 2)/std::pow(sigma_MKP1, 3);
            break;
        case 80:
            dJydsigma[80] = 1.0/sigma_Fos - 1.0*std::pow(Fos - mFos, 2)/std::pow(sigma_Fos, 3);
            break;
        case 81:
            dJydsigma[81] = 1.0/sigma_Jun - 1.0*std::pow(Jun - mJun, 2)/std::pow(sigma_Jun, 3);
            break;
        case 82:
            dJydsigma[82] = 1.0/sigma_Myc - 1.0*std::pow(Myc - mMyc, 2)/std::pow(sigma_Myc, 3);
            break;
        case 83:
            dJydsigma[83] = 1.0/sigma_bCATENIN - 1.0*std::pow(bCATENIN - mbCATENIN, 2)/std::pow(sigma_bCATENIN, 3);
            break;
        case 84:
            dJydsigma[84] = 1.0/sigma_PTEN - 1.0*std::pow(PTEN - mPTEN, 2)/std::pow(sigma_PTEN, 3);
            break;
        case 85:
            dJydsigma[85] = 1.0/sigma_AKT - 1.0*std::pow(AKT - mAKT, 2)/std::pow(sigma_AKT, 3);
            break;
        case 86:
            dJydsigma[86] = 1.0/sigma_PDK1 - 1.0*std::pow(PDK1 - mPDK1, 2)/std::pow(sigma_PDK1, 3);
            break;
        case 87:
            dJydsigma[87] = 1.0/sigma_Rictor - 1.0*std::pow(Rictor - mRictor, 2)/std::pow(sigma_Rictor, 3);
            break;
        case 88:
            dJydsigma[88] = 1.0/sigma_mTOR - 1.0*std::pow(mTOR - mmTOR, 2)/std::pow(sigma_mTOR, 3);
            break;
        case 89:
            dJydsigma[89] = 1.0/sigma_GSK3b - 1.0*std::pow(GSK3b - mGSK3b, 2)/std::pow(sigma_GSK3b, 3);
            break;
        case 90:
            dJydsigma[90] = 1.0/sigma_TSC1 - 1.0*std::pow(TSC1 - mTSC1, 2)/std::pow(sigma_TSC1, 3);
            break;
        case 91:
            dJydsigma[91] = 1.0/sigma_TSC2 - 1.0*std::pow(TSC2 - mTSC2, 2)/std::pow(sigma_TSC2, 3);
            break;
        case 92:
            dJydsigma[92] = 1.0/sigma_PKC - 1.0*std::pow(PKC - mPKC, 2)/std::pow(sigma_PKC, 3);
            break;
        case 93:
            dJydsigma[93] = 1.0/sigma_RKIP - 1.0*std::pow(RKIP - mRKIP, 2)/std::pow(sigma_RKIP, 3);
            break;
        case 94:
            dJydsigma[94] = 1.0/sigma_ERK - 1.0*std::pow(ERK - mERK, 2)/std::pow(sigma_ERK, 3);
            break;
        case 95:
            dJydsigma[95] = 1.0/sigma_FOXO - 1.0*std::pow(FOXO - mFOXO, 2)/std::pow(sigma_FOXO, 3);
            break;
        case 96:
            dJydsigma[96] = 1.0/sigma_Rheb - 1.0*std::pow(Rheb - mRheb, 2)/std::pow(sigma_Rheb, 3);
            break;
        case 97:
            dJydsigma[97] = 1.0/sigma_Raptor - 1.0*std::pow(Raptor - mRaptor, 2)/std::pow(sigma_Raptor, 3);
            break;
        case 98:
            dJydsigma[98] = 1.0/sigma_S6K - 1.0*std::pow(S6K - mS6K, 2)/std::pow(sigma_S6K, 3);
            break;
        case 99:
            dJydsigma[99] = 1.0/sigma_EIF4EBP1 - 1.0*std::pow(EIF4EBP1 - mEIF4EBP1, 2)/std::pow(sigma_EIF4EBP1, 3);
            break;
        case 100:
            dJydsigma[100] = 1.0/sigma_SOS - 1.0*std::pow(SOS - mSOS, 2)/std::pow(sigma_SOS, 3);
            break;
        case 101:
            dJydsigma[101] = 1.0/sigma_EIF4E - 1.0*std::pow(EIF4E - mEIF4E, 2)/std::pow(sigma_EIF4E, 3);
            break;
    }
}

} // namespace model_SPARCED_u87i
} // namespace amici
