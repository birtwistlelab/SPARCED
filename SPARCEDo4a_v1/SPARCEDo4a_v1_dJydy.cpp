#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "SPARCEDo4a_v1_k.h"
#include "SPARCEDo4a_v1_y.h"
#include "SPARCEDo4a_v1_sigmay.h"
#include "SPARCEDo4a_v1_my.h"
#include "SPARCEDo4a_v1_dJydy.h"

namespace amici {
namespace model_SPARCEDo4a_v1 {

void dJydy_SPARCEDo4a_v1(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = (-1.0*mp53 + 1.0*p53)/std::pow(sigma_p53, 2);
            break;
        case 1:
            dJydy[0] = (1.0*Mdm2 - 1.0*mMdm2)/std::pow(sigma_Mdm2, 2);
            break;
        case 2:
            dJydy[0] = (1.0*Wip1 - 1.0*mWip1)/std::pow(sigma_Wip1, 2);
            break;
        case 3:
            dJydy[0] = (1.0*BRCA2 - 1.0*mBRCA2)/std::pow(sigma_BRCA2, 2);
            break;
        case 4:
            dJydy[0] = (1.0*MSH6 - 1.0*mMSH6)/std::pow(sigma_MSH6, 2);
            break;
        case 5:
            dJydy[0] = (1.0*MGMT - 1.0*mMGMT)/std::pow(sigma_MGMT, 2);
            break;
        case 6:
            dJydy[0] = (1.0*ARF - 1.0*mARF)/std::pow(sigma_ARF, 2);
            break;
        case 7:
            dJydy[0] = (1.0*MDM4 - 1.0*mMDM4)/std::pow(sigma_MDM4, 2);
            break;
        case 8:
            dJydy[0] = (1.0*ATM - 1.0*mATM)/std::pow(sigma_ATM, 2);
            break;
        case 9:
            dJydy[0] = (1.0*ATR - 1.0*mATR)/std::pow(sigma_ATR, 2);
            break;
        case 10:
            dJydy[0] = (1.0*RB - 1.0*mRB)/std::pow(sigma_RB, 2);
            break;
        case 11:
            dJydy[0] = (1.0*E2F - 1.0*mE2F)/std::pow(sigma_E2F, 2);
            break;
        case 12:
            dJydy[0] = (1.0*Cd - 1.0*mCd)/std::pow(sigma_Cd, 2);
            break;
        case 13:
            dJydy[0] = (1.0*Ce - 1.0*mCe)/std::pow(sigma_Ce, 2);
            break;
        case 14:
            dJydy[0] = (1.0*Skp2 - 1.0*mSkp2)/std::pow(sigma_Skp2, 2);
            break;
        case 15:
            dJydy[0] = (1.0*Pai - 1.0*mPai)/std::pow(sigma_Pai, 2);
            break;
        case 16:
            dJydy[0] = (1.0*Pei - 1.0*mPei)/std::pow(sigma_Pei, 2);
            break;
        case 17:
            dJydy[0] = (1.0*Pbi - 1.0*mPbi)/std::pow(sigma_Pbi, 2);
            break;
        case 18:
            dJydy[0] = (1.0*Ca - 1.0*mCa)/std::pow(sigma_Ca, 2);
            break;
        case 19:
            dJydy[0] = (-1.0*mp27 + 1.0*p27)/std::pow(sigma_p27, 2);
            break;
        case 20:
            dJydy[0] = (1.0*Cdh1a - 1.0*mCdh1a)/std::pow(sigma_Cdh1a, 2);
            break;
        case 21:
            dJydy[0] = (1.0*Cb - 1.0*mCb)/std::pow(sigma_Cb, 2);
            break;
        case 22:
            dJydy[0] = (1.0*Cdc20 - 1.0*mCdc20)/std::pow(sigma_Cdc20, 2);
            break;
        case 23:
            dJydy[0] = (1.0*Wee1 - 1.0*mWee1)/std::pow(sigma_Wee1, 2);
            break;
        case 24:
            dJydy[0] = (1.0*Chk1 - 1.0*mChk1)/std::pow(sigma_Chk1, 2);
            break;
        case 25:
            dJydy[0] = (-1.0*mp21 + 1.0*p21)/std::pow(sigma_p21, 2);
            break;
        case 26:
            dJydy[0] = (1.0*L - 1.0*mL)/std::pow(sigma_L, 2);
            break;
        case 27:
            dJydy[0] = (1.0*R - 1.0*mR)/std::pow(sigma_R, 2);
            break;
        case 28:
            dJydy[0] = (1.0*flip - 1.0*mflip)/std::pow(sigma_flip, 2);
            break;
        case 29:
            dJydy[0] = (1.0*C8 - 1.0*mC8)/std::pow(sigma_C8, 2);
            break;
        case 30:
            dJydy[0] = (1.0*Bar - 1.0*mBar)/std::pow(sigma_Bar, 2);
            break;
        case 31:
            dJydy[0] = (1.0*C3 - 1.0*mC3)/std::pow(sigma_C3, 2);
            break;
        case 32:
            dJydy[0] = (1.0*C6 - 1.0*mC6)/std::pow(sigma_C6, 2);
            break;
        case 33:
            dJydy[0] = (1.0*XIAP - 1.0*mXIAP)/std::pow(sigma_XIAP, 2);
            break;
        case 34:
            dJydy[0] = (1.0*PARP - 1.0*mPARP)/std::pow(sigma_PARP, 2);
            break;
        case 35:
            dJydy[0] = (1.0*Bid - 1.0*mBid)/std::pow(sigma_Bid, 2);
            break;
        case 36:
            dJydy[0] = (1.0*Bcl2 - 1.0*mBcl2)/std::pow(sigma_Bcl2, 2);
            break;
        case 37:
            dJydy[0] = (1.0*Bax - 1.0*mBax)/std::pow(sigma_Bax, 2);
            break;
        case 38:
            dJydy[0] = (1.0*CytoC - 1.0*mCytoC)/std::pow(sigma_CytoC, 2);
            break;
        case 39:
            dJydy[0] = (1.0*Smac - 1.0*mSmac)/std::pow(sigma_Smac, 2);
            break;
        case 40:
            dJydy[0] = (1.0*Apaf - 1.0*mApaf)/std::pow(sigma_Apaf, 2);
            break;
        case 41:
            dJydy[0] = (1.0*C9 - 1.0*mC9)/std::pow(sigma_C9, 2);
            break;
        case 42:
            dJydy[0] = (1.0*BAD - 1.0*mBAD)/std::pow(sigma_BAD, 2);
            break;
        case 43:
            dJydy[0] = (1.0*PUMA - 1.0*mPUMA)/std::pow(sigma_PUMA, 2);
            break;
        case 44:
            dJydy[0] = (1.0*NOXA - 1.0*mNOXA)/std::pow(sigma_NOXA, 2);
            break;
        case 45:
            dJydy[0] = (1.0*BIM - 1.0*mBIM)/std::pow(sigma_BIM, 2);
            break;
        case 46:
            dJydy[0] = (1.0*E - 1.0*mE)/std::pow(sigma_E, 2);
            break;
        case 47:
            dJydy[0] = (1.0*H - 1.0*mH)/std::pow(sigma_H, 2);
            break;
        case 48:
            dJydy[0] = (1.0*HGF - 1.0*mHGF)/std::pow(sigma_HGF, 2);
            break;
        case 49:
            dJydy[0] = (1.0*P - 1.0*mP)/std::pow(sigma_P, 2);
            break;
        case 50:
            dJydy[0] = (1.0*F - 1.0*mF)/std::pow(sigma_F, 2);
            break;
        case 51:
            dJydy[0] = (1.0*I - 1.0*mI)/std::pow(sigma_I, 2);
            break;
        case 52:
            dJydy[0] = (1.0*IN - 1.0*mIN)/std::pow(sigma_IN, 2);
            break;
        case 53:
            dJydy[0] = (1.0*E1 - 1.0*mE1)/std::pow(sigma_E1, 2);
            break;
        case 54:
            dJydy[0] = (1.0*E2 - 1.0*mE2)/std::pow(sigma_E2, 2);
            break;
        case 55:
            dJydy[0] = (1.0*E3 - 1.0*mE3)/std::pow(sigma_E3, 2);
            break;
        case 56:
            dJydy[0] = (1.0*E4 - 1.0*mE4)/std::pow(sigma_E4, 2);
            break;
        case 57:
            dJydy[0] = (1.0*Ev3 - 1.0*mEv3)/std::pow(sigma_Ev3, 2);
            break;
        case 58:
            dJydy[0] = (1.0*Met - 1.0*mMet)/std::pow(sigma_Met, 2);
            break;
        case 59:
            dJydy[0] = (1.0*Pr - 1.0*mPr)/std::pow(sigma_Pr, 2);
            break;
        case 60:
            dJydy[0] = (1.0*Fr - 1.0*mFr)/std::pow(sigma_Fr, 2);
            break;
        case 61:
            dJydy[0] = (1.0*Ir - 1.0*mIr)/std::pow(sigma_Ir, 2);
            break;
        case 62:
            dJydy[0] = (1.0*Isr - 1.0*mIsr)/std::pow(sigma_Isr, 2);
            break;
        case 63:
            dJydy[0] = (1.0*IRS - 1.0*mIRS)/std::pow(sigma_IRS, 2);
            break;
        case 64:
            dJydy[0] = (1.0*Sp - 1.0*mSp)/std::pow(sigma_Sp, 2);
            break;
        case 65:
            dJydy[0] = (1.0*Cbl - 1.0*mCbl)/std::pow(sigma_Cbl, 2);
            break;
        case 66:
            dJydy[0] = (1.0*G2 - 1.0*mG2)/std::pow(sigma_G2, 2);
            break;
        case 67:
            dJydy[0] = (1.0*PLCg - 1.0*mPLCg)/std::pow(sigma_PLCg, 2);
            break;
        case 68:
            dJydy[0] = (1.0*PI3KC1 - 1.0*mPI3KC1)/std::pow(sigma_PI3KC1, 2);
            break;
        case 69:
            dJydy[0] = (1.0*PI3KR1 - 1.0*mPI3KR1)/std::pow(sigma_PI3KR1, 2);
            break;
        case 70:
            dJydy[0] = (1.0*PI3K2 - 1.0*mPI3K2)/std::pow(sigma_PI3K2, 2);
            break;
        case 71:
            dJydy[0] = (1.0*GRP - 1.0*mGRP)/std::pow(sigma_GRP, 2);
            break;
        case 72:
            dJydy[0] = (1.0*Ras - 1.0*mRas)/std::pow(sigma_Ras, 2);
            break;
        case 73:
            dJydy[0] = (1.0*NF1 - 1.0*mNF1)/std::pow(sigma_NF1, 2);
            break;
        case 74:
            dJydy[0] = (1.0*CRaf - 1.0*mCRaf)/std::pow(sigma_CRaf, 2);
            break;
        case 75:
            dJydy[0] = (1.0*BRaf - 1.0*mBRaf)/std::pow(sigma_BRaf, 2);
            break;
        case 76:
            dJydy[0] = (1.0*MEK - 1.0*mMEK)/std::pow(sigma_MEK, 2);
            break;
        case 77:
            dJydy[0] = (1.0*MKP3 - 1.0*mMKP3)/std::pow(sigma_MKP3, 2);
            break;
        case 78:
            dJydy[0] = (1.0*RSK - 1.0*mRSK)/std::pow(sigma_RSK, 2);
            break;
        case 79:
            dJydy[0] = (1.0*MKP1 - 1.0*mMKP1)/std::pow(sigma_MKP1, 2);
            break;
        case 80:
            dJydy[0] = (1.0*Fos - 1.0*mFos)/std::pow(sigma_Fos, 2);
            break;
        case 81:
            dJydy[0] = (1.0*Jun - 1.0*mJun)/std::pow(sigma_Jun, 2);
            break;
        case 82:
            dJydy[0] = (1.0*Myc - 1.0*mMyc)/std::pow(sigma_Myc, 2);
            break;
        case 83:
            dJydy[0] = (1.0*bCATENIN - 1.0*mbCATENIN)/std::pow(sigma_bCATENIN, 2);
            break;
        case 84:
            dJydy[0] = (1.0*PTEN - 1.0*mPTEN)/std::pow(sigma_PTEN, 2);
            break;
        case 85:
            dJydy[0] = (1.0*AKT - 1.0*mAKT)/std::pow(sigma_AKT, 2);
            break;
        case 86:
            dJydy[0] = (1.0*PDK1 - 1.0*mPDK1)/std::pow(sigma_PDK1, 2);
            break;
        case 87:
            dJydy[0] = (1.0*Rictor - 1.0*mRictor)/std::pow(sigma_Rictor, 2);
            break;
        case 88:
            dJydy[0] = (1.0*mTOR - 1.0*mmTOR)/std::pow(sigma_mTOR, 2);
            break;
        case 89:
            dJydy[0] = (1.0*GSK3b - 1.0*mGSK3b)/std::pow(sigma_GSK3b, 2);
            break;
        case 90:
            dJydy[0] = (1.0*TSC1 - 1.0*mTSC1)/std::pow(sigma_TSC1, 2);
            break;
        case 91:
            dJydy[0] = (1.0*TSC2 - 1.0*mTSC2)/std::pow(sigma_TSC2, 2);
            break;
        case 92:
            dJydy[0] = (1.0*PKC - 1.0*mPKC)/std::pow(sigma_PKC, 2);
            break;
        case 93:
            dJydy[0] = (1.0*RKIP - 1.0*mRKIP)/std::pow(sigma_RKIP, 2);
            break;
        case 94:
            dJydy[0] = (1.0*ERK - 1.0*mERK)/std::pow(sigma_ERK, 2);
            break;
        case 95:
            dJydy[0] = (1.0*FOXO - 1.0*mFOXO)/std::pow(sigma_FOXO, 2);
            break;
        case 96:
            dJydy[0] = (1.0*Rheb - 1.0*mRheb)/std::pow(sigma_Rheb, 2);
            break;
        case 97:
            dJydy[0] = (1.0*Raptor - 1.0*mRaptor)/std::pow(sigma_Raptor, 2);
            break;
        case 98:
            dJydy[0] = (1.0*S6K - 1.0*mS6K)/std::pow(sigma_S6K, 2);
            break;
        case 99:
            dJydy[0] = (1.0*EIF4EBP1 - 1.0*mEIF4EBP1)/std::pow(sigma_EIF4EBP1, 2);
            break;
        case 100:
            dJydy[0] = (1.0*SOS - 1.0*mSOS)/std::pow(sigma_SOS, 2);
            break;
        case 101:
            dJydy[0] = (1.0*EIF4E - 1.0*mEIF4E)/std::pow(sigma_EIF4E, 2);
            break;
    }
}

} // namespace model_SPARCEDo4a_v1
} // namespace amici
