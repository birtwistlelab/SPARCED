#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <array>

#include "y.h"
#include "my.h"
#include "p.h"
#include "k.h"
#include "sigmay.h"
#include "dJydy.h"

namespace amici {
namespace model_SPARCED_tutorial {

void dJydy_SPARCED_tutorial(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = (-1.0*myRibosome + 1.0*yRibosome)/std::pow(sigma_yRibosome, 2);
            break;
        case 1:
            dJydy[0] = (-1.0*myp53inac + 1.0*yp53inac)/std::pow(sigma_yp53inac, 2);
            break;
        case 2:
            dJydy[0] = (-1.0*myp53ac + 1.0*yp53ac)/std::pow(sigma_yp53ac, 2);
            break;
        case 3:
            dJydy[0] = (-1.0*myMDM2 + 1.0*yMDM2)/std::pow(sigma_yMDM2, 2);
            break;
        case 4:
            dJydy[0] = (-1.0*myWip1 + 1.0*yWip1)/std::pow(sigma_yWip1, 2);
            break;
        case 5:
            dJydy[0] = (-1.0*myATMP + 1.0*yATMP)/std::pow(sigma_yATMP, 2);
            break;
        case 6:
            dJydy[0] = (-1.0*myATRac + 1.0*yATRac)/std::pow(sigma_yATRac, 2);
            break;
        case 7:
            dJydy[0] = (-1.0*myMDM2product1 + 1.0*yMDM2product1)/std::pow(sigma_yMDM2product1, 2);
            break;
        case 8:
            dJydy[0] = (-1.0*myMDM2product2 + 1.0*yMDM2product2)/std::pow(sigma_yMDM2product2, 2);
            break;
        case 9:
            dJydy[0] = (-1.0*myMDM2product3 + 1.0*yMDM2product3)/std::pow(sigma_yMDM2product3, 2);
            break;
        case 10:
            dJydy[0] = (-1.0*myMDM2product4 + 1.0*yMDM2product4)/std::pow(sigma_yMDM2product4, 2);
            break;
        case 11:
            dJydy[0] = (-1.0*myMDM2product5 + 1.0*yMDM2product5)/std::pow(sigma_yMDM2product5, 2);
            break;
        case 12:
            dJydy[0] = (-1.0*myMDM2product6 + 1.0*yMDM2product6)/std::pow(sigma_yMDM2product6, 2);
            break;
        case 13:
            dJydy[0] = (-1.0*myMDM2product7 + 1.0*yMDM2product7)/std::pow(sigma_yMDM2product7, 2);
            break;
        case 14:
            dJydy[0] = (-1.0*myMDM2product8 + 1.0*yMDM2product8)/std::pow(sigma_yMDM2product8, 2);
            break;
        case 15:
            dJydy[0] = (-1.0*myMDM2product9 + 1.0*yMDM2product9)/std::pow(sigma_yMDM2product9, 2);
            break;
        case 16:
            dJydy[0] = (-1.0*myMDM2pro + 1.0*yMDM2pro)/std::pow(sigma_yMDM2pro, 2);
            break;
        case 17:
            dJydy[0] = (-1.0*myWip1product1 + 1.0*yWip1product1)/std::pow(sigma_yWip1product1, 2);
            break;
        case 18:
            dJydy[0] = (-1.0*myWip1product2 + 1.0*yWip1product2)/std::pow(sigma_yWip1product2, 2);
            break;
        case 19:
            dJydy[0] = (-1.0*myWip1product3 + 1.0*yWip1product3)/std::pow(sigma_yWip1product3, 2);
            break;
        case 20:
            dJydy[0] = (-1.0*myWip1product4 + 1.0*yWip1product4)/std::pow(sigma_yWip1product4, 2);
            break;
        case 21:
            dJydy[0] = (-1.0*myWip1product5 + 1.0*yWip1product5)/std::pow(sigma_yWip1product5, 2);
            break;
        case 22:
            dJydy[0] = (-1.0*myWip1product6 + 1.0*yWip1product6)/std::pow(sigma_yWip1product6, 2);
            break;
        case 23:
            dJydy[0] = (-1.0*myWip1product7 + 1.0*yWip1product7)/std::pow(sigma_yWip1product7, 2);
            break;
        case 24:
            dJydy[0] = (-1.0*myWip1product8 + 1.0*yWip1product8)/std::pow(sigma_yWip1product8, 2);
            break;
        case 25:
            dJydy[0] = (-1.0*myWip1product9 + 1.0*yWip1product9)/std::pow(sigma_yWip1product9, 2);
            break;
        case 26:
            dJydy[0] = (-1.0*myWip1pro + 1.0*yWip1pro)/std::pow(sigma_yWip1pro, 2);
            break;
        case 27:
            dJydy[0] = (-1.0*myBRCA2 + 1.0*yBRCA2)/std::pow(sigma_yBRCA2, 2);
            break;
        case 28:
            dJydy[0] = (-1.0*myMSH6 + 1.0*yMSH6)/std::pow(sigma_yMSH6, 2);
            break;
        case 29:
            dJydy[0] = (-1.0*myMGMT + 1.0*yMGMT)/std::pow(sigma_yMGMT, 2);
            break;
        case 30:
            dJydy[0] = (-1.0*mydamageDSB + 1.0*ydamageDSB)/std::pow(sigma_ydamageDSB, 2);
            break;
        case 31:
            dJydy[0] = (-1.0*mydamageSSB + 1.0*ydamageSSB)/std::pow(sigma_ydamageSSB, 2);
            break;
        case 32:
            dJydy[0] = (-1.0*myppAKT_MDM2 + 1.0*yppAKT_MDM2)/std::pow(sigma_yppAKT_MDM2, 2);
            break;
        case 33:
            dJydy[0] = (-1.0*mypMDM2 + 1.0*ypMDM2)/std::pow(sigma_ypMDM2, 2);
            break;
        case 34:
            dJydy[0] = (-1.0*myARF + 1.0*yARF)/std::pow(sigma_yARF, 2);
            break;
        case 35:
            dJydy[0] = (-1.0*myMDM4 + 1.0*yMDM4)/std::pow(sigma_yMDM4, 2);
            break;
        case 36:
            dJydy[0] = (-1.0*myp53ac_MDM4 + 1.0*yp53ac_MDM4)/std::pow(sigma_yp53ac_MDM4, 2);
            break;
        case 37:
            dJydy[0] = (-1.0*myATMinac + 1.0*yATMinac)/std::pow(sigma_yATMinac, 2);
            break;
        case 38:
            dJydy[0] = (-1.0*myATRinac + 1.0*yATRinac)/std::pow(sigma_yATRinac, 2);
            break;
        case 39:
            dJydy[0] = (-1.0*mypRB + 1.0*ypRB)/std::pow(sigma_ypRB, 2);
            break;
        case 40:
            dJydy[0] = (-1.0*mypRBp + 1.0*ypRBp)/std::pow(sigma_ypRBp, 2);
            break;
        case 41:
            dJydy[0] = (-1.0*mypRBpp + 1.0*ypRBpp)/std::pow(sigma_ypRBpp, 2);
            break;
        case 42:
            dJydy[0] = (-1.0*myE2F + 1.0*yE2F)/std::pow(sigma_yE2F, 2);
            break;
        case 43:
            dJydy[0] = (-1.0*myCd + 1.0*yCd)/std::pow(sigma_yCd, 2);
            break;
        case 44:
            dJydy[0] = (-1.0*myMdi + 1.0*yMdi)/std::pow(sigma_yMdi, 2);
            break;
        case 45:
            dJydy[0] = (-1.0*myMd + 1.0*yMd)/std::pow(sigma_yMd, 2);
            break;
        case 46:
            dJydy[0] = (-1.0*myMdp27 + 1.0*yMdp27)/std::pow(sigma_yMdp27, 2);
            break;
        case 47:
            dJydy[0] = (-1.0*myCe + 1.0*yCe)/std::pow(sigma_yCe, 2);
            break;
        case 48:
            dJydy[0] = (-1.0*myMei + 1.0*yMei)/std::pow(sigma_yMei, 2);
            break;
        case 49:
            dJydy[0] = (-1.0*myMe + 1.0*yMe)/std::pow(sigma_yMe, 2);
            break;
        case 50:
            dJydy[0] = (-1.0*mySkp2 + 1.0*ySkp2)/std::pow(sigma_ySkp2, 2);
            break;
        case 51:
            dJydy[0] = (-1.0*myMep27 + 1.0*yMep27)/std::pow(sigma_yMep27, 2);
            break;
        case 52:
            dJydy[0] = (-1.0*myPe + 1.0*yPe)/std::pow(sigma_yPe, 2);
            break;
        case 53:
            dJydy[0] = (-1.0*myPai + 1.0*yPai)/std::pow(sigma_yPai, 2);
            break;
        case 54:
            dJydy[0] = (-1.0*myPei + 1.0*yPei)/std::pow(sigma_yPei, 2);
            break;
        case 55:
            dJydy[0] = (-1.0*myPbi + 1.0*yPbi)/std::pow(sigma_yPbi, 2);
            break;
        case 56:
            dJydy[0] = (-1.0*myCa + 1.0*yCa)/std::pow(sigma_yCa, 2);
            break;
        case 57:
            dJydy[0] = (-1.0*myMai + 1.0*yMai)/std::pow(sigma_yMai, 2);
            break;
        case 58:
            dJydy[0] = (-1.0*myMa + 1.0*yMa)/std::pow(sigma_yMa, 2);
            break;
        case 59:
            dJydy[0] = (-1.0*myMap27 + 1.0*yMap27)/std::pow(sigma_yMap27, 2);
            break;
        case 60:
            dJydy[0] = (-1.0*myp27 + 1.0*yp27)/std::pow(sigma_yp27, 2);
            break;
        case 61:
            dJydy[0] = (-1.0*myCdh1i + 1.0*yCdh1i)/std::pow(sigma_yCdh1i, 2);
            break;
        case 62:
            dJydy[0] = (-1.0*myCdh1a + 1.0*yCdh1a)/std::pow(sigma_yCdh1a, 2);
            break;
        case 63:
            dJydy[0] = (-1.0*myE2Fp + 1.0*yE2Fp)/std::pow(sigma_yE2Fp, 2);
            break;
        case 64:
            dJydy[0] = (-1.0*myp27p + 1.0*yp27p)/std::pow(sigma_yp27p, 2);
            break;
        case 65:
            dJydy[0] = (-1.0*myPa + 1.0*yPa)/std::pow(sigma_yPa, 2);
            break;
        case 66:
            dJydy[0] = (-1.0*myCb + 1.0*yCb)/std::pow(sigma_yCb, 2);
            break;
        case 67:
            dJydy[0] = (-1.0*myMbi + 1.0*yMbi)/std::pow(sigma_yMbi, 2);
            break;
        case 68:
            dJydy[0] = (-1.0*myMb + 1.0*yMb)/std::pow(sigma_yMb, 2);
            break;
        case 69:
            dJydy[0] = (-1.0*myCdc20i + 1.0*yCdc20i)/std::pow(sigma_yCdc20i, 2);
            break;
        case 70:
            dJydy[0] = (-1.0*myCdc20a + 1.0*yCdc20a)/std::pow(sigma_yCdc20a, 2);
            break;
        case 71:
            dJydy[0] = (-1.0*myPb + 1.0*yPb)/std::pow(sigma_yPb, 2);
            break;
        case 72:
            dJydy[0] = (-1.0*myWee1 + 1.0*yWee1)/std::pow(sigma_yWee1, 2);
            break;
        case 73:
            dJydy[0] = (-1.0*myWee1p + 1.0*yWee1p)/std::pow(sigma_yWee1p, 2);
            break;
        case 74:
            dJydy[0] = (-1.0*myMbp27 + 1.0*yMbp27)/std::pow(sigma_yMbp27, 2);
            break;
        case 75:
            dJydy[0] = (-1.0*myChk1 + 1.0*yChk1)/std::pow(sigma_yChk1, 2);
            break;
        case 76:
            dJydy[0] = (-1.0*mypRBc1 + 1.0*ypRBc1)/std::pow(sigma_ypRBc1, 2);
            break;
        case 77:
            dJydy[0] = (-1.0*mypRBc2 + 1.0*ypRBc2)/std::pow(sigma_ypRBc2, 2);
            break;
        case 78:
            dJydy[0] = (-1.0*myp21 + 1.0*yp21)/std::pow(sigma_yp21, 2);
            break;
        case 79:
            dJydy[0] = (-1.0*myMdp21 + 1.0*yMdp21)/std::pow(sigma_yMdp21, 2);
            break;
        case 80:
            dJydy[0] = (-1.0*myMep21 + 1.0*yMep21)/std::pow(sigma_yMep21, 2);
            break;
        case 81:
            dJydy[0] = (-1.0*myMap21 + 1.0*yMap21)/std::pow(sigma_yMap21, 2);
            break;
        case 82:
            dJydy[0] = (-1.0*myMbp21 + 1.0*yMbp21)/std::pow(sigma_yMbp21, 2);
            break;
        case 83:
            dJydy[0] = (-1.0*myL + 1.0*yL)/std::pow(sigma_yL, 2);
            break;
        case 84:
            dJydy[0] = (-1.0*myR + 1.0*yR)/std::pow(sigma_yR, 2);
            break;
        case 85:
            dJydy[0] = (-1.0*myL_R + 1.0*yL_R)/std::pow(sigma_yL_R, 2);
            break;
        case 86:
            dJydy[0] = (-1.0*myRactive + 1.0*yRactive)/std::pow(sigma_yRactive, 2);
            break;
        case 87:
            dJydy[0] = (-1.0*myflip + 1.0*yflip)/std::pow(sigma_yflip, 2);
            break;
        case 88:
            dJydy[0] = (-1.0*myRactive_flip + 1.0*yRactive_flip)/std::pow(sigma_yRactive_flip, 2);
            break;
        case 89:
            dJydy[0] = (-1.0*mypC8 + 1.0*ypC8)/std::pow(sigma_ypC8, 2);
            break;
        case 90:
            dJydy[0] = (-1.0*myRactive_pC8 + 1.0*yRactive_pC8)/std::pow(sigma_yRactive_pC8, 2);
            break;
        case 91:
            dJydy[0] = (-1.0*myC8 + 1.0*yC8)/std::pow(sigma_yC8, 2);
            break;
        case 92:
            dJydy[0] = (-1.0*myBar + 1.0*yBar)/std::pow(sigma_yBar, 2);
            break;
        case 93:
            dJydy[0] = (-1.0*myC8_Bar + 1.0*yC8_Bar)/std::pow(sigma_yC8_Bar, 2);
            break;
        case 94:
            dJydy[0] = (-1.0*mypC3 + 1.0*ypC3)/std::pow(sigma_ypC3, 2);
            break;
        case 95:
            dJydy[0] = (-1.0*myC8_pC3 + 1.0*yC8_pC3)/std::pow(sigma_yC8_pC3, 2);
            break;
        case 96:
            dJydy[0] = (-1.0*myC3 + 1.0*yC3)/std::pow(sigma_yC3, 2);
            break;
        case 97:
            dJydy[0] = (-1.0*mypC6 + 1.0*ypC6)/std::pow(sigma_ypC6, 2);
            break;
        case 98:
            dJydy[0] = (-1.0*myC3_pC6 + 1.0*yC3_pC6)/std::pow(sigma_yC3_pC6, 2);
            break;
        case 99:
            dJydy[0] = (-1.0*myC6 + 1.0*yC6)/std::pow(sigma_yC6, 2);
            break;
        case 100:
            dJydy[0] = (-1.0*myC6_C8 + 1.0*yC6_C8)/std::pow(sigma_yC6_C8, 2);
            break;
        case 101:
            dJydy[0] = (-1.0*myXIAP + 1.0*yXIAP)/std::pow(sigma_yXIAP, 2);
            break;
        case 102:
            dJydy[0] = (-1.0*myC3_XIAP + 1.0*yC3_XIAP)/std::pow(sigma_yC3_XIAP, 2);
            break;
        case 103:
            dJydy[0] = (-1.0*myPARP + 1.0*yPARP)/std::pow(sigma_yPARP, 2);
            break;
        case 104:
            dJydy[0] = (-1.0*myC3_PARP + 1.0*yC3_PARP)/std::pow(sigma_yC3_PARP, 2);
            break;
        case 105:
            dJydy[0] = (-1.0*mycPARP + 1.0*ycPARP)/std::pow(sigma_ycPARP, 2);
            break;
        case 106:
            dJydy[0] = (-1.0*myBid + 1.0*yBid)/std::pow(sigma_yBid, 2);
            break;
        case 107:
            dJydy[0] = (-1.0*myC8_Bid + 1.0*yC8_Bid)/std::pow(sigma_yC8_Bid, 2);
            break;
        case 108:
            dJydy[0] = (-1.0*mytBid + 1.0*ytBid)/std::pow(sigma_ytBid, 2);
            break;
        case 109:
            dJydy[0] = (-1.0*myBcl2c + 1.0*yBcl2c)/std::pow(sigma_yBcl2c, 2);
            break;
        case 110:
            dJydy[0] = (-1.0*mytBid_Bcl2c + 1.0*ytBid_Bcl2c)/std::pow(sigma_ytBid_Bcl2c, 2);
            break;
        case 111:
            dJydy[0] = (-1.0*myBax + 1.0*yBax)/std::pow(sigma_yBax, 2);
            break;
        case 112:
            dJydy[0] = (-1.0*mytBid_Bax + 1.0*ytBid_Bax)/std::pow(sigma_ytBid_Bax, 2);
            break;
        case 113:
            dJydy[0] = (-1.0*myBaxactive + 1.0*yBaxactive)/std::pow(sigma_yBaxactive, 2);
            break;
        case 114:
            dJydy[0] = (-1.0*myBaxm + 1.0*yBaxm)/std::pow(sigma_yBaxm, 2);
            break;
        case 115:
            dJydy[0] = (-1.0*myBcl2 + 1.0*yBcl2)/std::pow(sigma_yBcl2, 2);
            break;
        case 116:
            dJydy[0] = (-1.0*myBaxm_Bcl2 + 1.0*yBaxm_Bcl2)/std::pow(sigma_yBaxm_Bcl2, 2);
            break;
        case 117:
            dJydy[0] = (-1.0*myBax2 + 1.0*yBax2)/std::pow(sigma_yBax2, 2);
            break;
        case 118:
            dJydy[0] = (-1.0*myBax2_Bcl2 + 1.0*yBax2_Bcl2)/std::pow(sigma_yBax2_Bcl2, 2);
            break;
        case 119:
            dJydy[0] = (-1.0*myBax4 + 1.0*yBax4)/std::pow(sigma_yBax4, 2);
            break;
        case 120:
            dJydy[0] = (-1.0*myBax4_Bcl2 + 1.0*yBax4_Bcl2)/std::pow(sigma_yBax4_Bcl2, 2);
            break;
        case 121:
            dJydy[0] = (-1.0*myM + 1.0*yM)/std::pow(sigma_yM, 2);
            break;
        case 122:
            dJydy[0] = (-1.0*myBax4_M + 1.0*yBax4_M)/std::pow(sigma_yBax4_M, 2);
            break;
        case 123:
            dJydy[0] = (-1.0*myMactive + 1.0*yMactive)/std::pow(sigma_yMactive, 2);
            break;
        case 124:
            dJydy[0] = (-1.0*myCytoCm + 1.0*yCytoCm)/std::pow(sigma_yCytoCm, 2);
            break;
        case 125:
            dJydy[0] = (-1.0*myMactive_CytoCm + 1.0*yMactive_CytoCm)/std::pow(sigma_yMactive_CytoCm, 2);
            break;
        case 126:
            dJydy[0] = (-1.0*myCytoCr + 1.0*yCytoCr)/std::pow(sigma_yCytoCr, 2);
            break;
        case 127:
            dJydy[0] = (-1.0*mySmacm + 1.0*ySmacm)/std::pow(sigma_ySmacm, 2);
            break;
        case 128:
            dJydy[0] = (-1.0*myMactive_Smacm + 1.0*yMactive_Smacm)/std::pow(sigma_yMactive_Smacm, 2);
            break;
        case 129:
            dJydy[0] = (-1.0*mySmacr + 1.0*ySmacr)/std::pow(sigma_ySmacr, 2);
            break;
        case 130:
            dJydy[0] = (-1.0*myCytoC + 1.0*yCytoC)/std::pow(sigma_yCytoC, 2);
            break;
        case 131:
            dJydy[0] = (-1.0*myApaf + 1.0*yApaf)/std::pow(sigma_yApaf, 2);
            break;
        case 132:
            dJydy[0] = (-1.0*myCytoC_Apaf + 1.0*yCytoC_Apaf)/std::pow(sigma_yCytoC_Apaf, 2);
            break;
        case 133:
            dJydy[0] = (-1.0*myApafactive + 1.0*yApafactive)/std::pow(sigma_yApafactive, 2);
            break;
        case 134:
            dJydy[0] = (-1.0*mypC9 + 1.0*ypC9)/std::pow(sigma_ypC9, 2);
            break;
        case 135:
            dJydy[0] = (-1.0*myApop + 1.0*yApop)/std::pow(sigma_yApop, 2);
            break;
        case 136:
            dJydy[0] = (-1.0*myApop_C3 + 1.0*yApop_C3)/std::pow(sigma_yApop_C3, 2);
            break;
        case 137:
            dJydy[0] = (-1.0*mySmac + 1.0*ySmac)/std::pow(sigma_ySmac, 2);
            break;
        case 138:
            dJydy[0] = (-1.0*myApop_XIAP + 1.0*yApop_XIAP)/std::pow(sigma_yApop_XIAP, 2);
            break;
        case 139:
            dJydy[0] = (-1.0*mySmac_XIAP + 1.0*ySmac_XIAP)/std::pow(sigma_ySmac_XIAP, 2);
            break;
        case 140:
            dJydy[0] = (-1.0*myC3_Ub + 1.0*yC3_Ub)/std::pow(sigma_yC3_Ub, 2);
            break;
        case 141:
            dJydy[0] = (-1.0*myBAD + 1.0*yBAD)/std::pow(sigma_yBAD, 2);
            break;
        case 142:
            dJydy[0] = (-1.0*myPUMA + 1.0*yPUMA)/std::pow(sigma_yPUMA, 2);
            break;
        case 143:
            dJydy[0] = (-1.0*myNOXA + 1.0*yNOXA)/std::pow(sigma_yNOXA, 2);
            break;
        case 144:
            dJydy[0] = (-1.0*myBcl2c_BAD + 1.0*yBcl2c_BAD)/std::pow(sigma_yBcl2c_BAD, 2);
            break;
        case 145:
            dJydy[0] = (-1.0*myBcl2c_PUMA + 1.0*yBcl2c_PUMA)/std::pow(sigma_yBcl2c_PUMA, 2);
            break;
        case 146:
            dJydy[0] = (-1.0*myBcl2c_NOXA + 1.0*yBcl2c_NOXA)/std::pow(sigma_yBcl2c_NOXA, 2);
            break;
        case 147:
            dJydy[0] = (-1.0*myBIM + 1.0*yBIM)/std::pow(sigma_yBIM, 2);
            break;
        case 148:
            dJydy[0] = (-1.0*myBIM_Bax + 1.0*yBIM_Bax)/std::pow(sigma_yBIM_Bax, 2);
            break;
        case 149:
            dJydy[0] = (-1.0*myBcl2c_BIM + 1.0*yBcl2c_BIM)/std::pow(sigma_yBcl2c_BIM, 2);
            break;
        case 150:
            dJydy[0] = (-1.0*myppERK_BIM + 1.0*yppERK_BIM)/std::pow(sigma_yppERK_BIM, 2);
            break;
        case 151:
            dJydy[0] = (-1.0*mypBIM + 1.0*ypBIM)/std::pow(sigma_ypBIM, 2);
            break;
        case 152:
            dJydy[0] = (-1.0*myppAKT_BAD + 1.0*yppAKT_BAD)/std::pow(sigma_yppAKT_BAD, 2);
            break;
        case 153:
            dJydy[0] = (-1.0*mypBAD + 1.0*ypBAD)/std::pow(sigma_ypBAD, 2);
            break;
        case 154:
            dJydy[0] = (-1.0*myppERK_BAD + 1.0*yppERK_BAD)/std::pow(sigma_yppERK_BAD, 2);
            break;
        case 155:
            dJydy[0] = (-1.0*myE + 1.0*yE)/std::pow(sigma_yE, 2);
            break;
        case 156:
            dJydy[0] = (-1.0*myH + 1.0*yH)/std::pow(sigma_yH, 2);
            break;
        case 157:
            dJydy[0] = (-1.0*myHGF + 1.0*yHGF)/std::pow(sigma_yHGF, 2);
            break;
        case 158:
            dJydy[0] = (-1.0*myP + 1.0*yP)/std::pow(sigma_yP, 2);
            break;
        case 159:
            dJydy[0] = (-1.0*myF + 1.0*yF)/std::pow(sigma_yF, 2);
            break;
        case 160:
            dJydy[0] = (-1.0*myI + 1.0*yI)/std::pow(sigma_yI, 2);
            break;
        case 161:
            dJydy[0] = (-1.0*myINS + 1.0*yINS)/std::pow(sigma_yINS, 2);
            break;
        case 162:
            dJydy[0] = (-1.0*myE1 + 1.0*yE1)/std::pow(sigma_yE1, 2);
            break;
        case 163:
            dJydy[0] = (-1.0*mypE1 + 1.0*ypE1)/std::pow(sigma_ypE1, 2);
            break;
        case 164:
            dJydy[0] = (-1.0*myE2 + 1.0*yE2)/std::pow(sigma_yE2, 2);
            break;
        case 165:
            dJydy[0] = (-1.0*mypE2 + 1.0*ypE2)/std::pow(sigma_ypE2, 2);
            break;
        case 166:
            dJydy[0] = (-1.0*myE3 + 1.0*yE3)/std::pow(sigma_yE3, 2);
            break;
        case 167:
            dJydy[0] = (-1.0*myE4 + 1.0*yE4)/std::pow(sigma_yE4, 2);
            break;
        case 168:
            dJydy[0] = (-1.0*mypE4 + 1.0*ypE4)/std::pow(sigma_ypE4, 2);
            break;
        case 169:
            dJydy[0] = (-1.0*myEv3 + 1.0*yEv3)/std::pow(sigma_yEv3, 2);
            break;
        case 170:
            dJydy[0] = (-1.0*myMet + 1.0*yMet)/std::pow(sigma_yMet, 2);
            break;
        case 171:
            dJydy[0] = (-1.0*myPr + 1.0*yPr)/std::pow(sigma_yPr, 2);
            break;
        case 172:
            dJydy[0] = (-1.0*myFr + 1.0*yFr)/std::pow(sigma_yFr, 2);
            break;
        case 173:
            dJydy[0] = (-1.0*myIr + 1.0*yIr)/std::pow(sigma_yIr, 2);
            break;
        case 174:
            dJydy[0] = (-1.0*myIsr + 1.0*yIsr)/std::pow(sigma_yIsr, 2);
            break;
        case 175:
            dJydy[0] = (-1.0*myE1E1 + 1.0*yE1E1)/std::pow(sigma_yE1E1, 2);
            break;
        case 176:
            dJydy[0] = (-1.0*myE1E2 + 1.0*yE1E2)/std::pow(sigma_yE1E2, 2);
            break;
        case 177:
            dJydy[0] = (-1.0*myE1E3 + 1.0*yE1E3)/std::pow(sigma_yE1E3, 2);
            break;
        case 178:
            dJydy[0] = (-1.0*myE1E4 + 1.0*yE1E4)/std::pow(sigma_yE1E4, 2);
            break;
        case 179:
            dJydy[0] = (-1.0*myE2E2 + 1.0*yE2E2)/std::pow(sigma_yE2E2, 2);
            break;
        case 180:
            dJydy[0] = (-1.0*myE2E3 + 1.0*yE2E3)/std::pow(sigma_yE2E3, 2);
            break;
        case 181:
            dJydy[0] = (-1.0*myE2E4 + 1.0*yE2E4)/std::pow(sigma_yE2E4, 2);
            break;
        case 182:
            dJydy[0] = (-1.0*myE3E4 + 1.0*yE3E4)/std::pow(sigma_yE3E4, 2);
            break;
        case 183:
            dJydy[0] = (-1.0*myE4E4 + 1.0*yE4E4)/std::pow(sigma_yE4E4, 2);
            break;
        case 184:
            dJydy[0] = (-1.0*myMet_Met + 1.0*yMet_Met)/std::pow(sigma_yMet_Met, 2);
            break;
        case 185:
            dJydy[0] = (-1.0*myFrFr + 1.0*yFrFr)/std::pow(sigma_yFrFr, 2);
            break;
        case 186:
            dJydy[0] = (-1.0*myIrIr + 1.0*yIrIr)/std::pow(sigma_yIrIr, 2);
            break;
        case 187:
            dJydy[0] = (-1.0*myIsr_Isr + 1.0*yIsr_Isr)/std::pow(sigma_yIsr_Isr, 2);
            break;
        case 188:
            dJydy[0] = (-1.0*myEE1 + 1.0*yEE1)/std::pow(sigma_yEE1, 2);
            break;
        case 189:
            dJydy[0] = (-1.0*myHE3 + 1.0*yHE3)/std::pow(sigma_yHE3, 2);
            break;
        case 190:
            dJydy[0] = (-1.0*myHE4 + 1.0*yHE4)/std::pow(sigma_yHE4, 2);
            break;
        case 191:
            dJydy[0] = (-1.0*myHGF_Met + 1.0*yHGF_Met)/std::pow(sigma_yHGF_Met, 2);
            break;
        case 192:
            dJydy[0] = (-1.0*myPPr + 1.0*yPPr)/std::pow(sigma_yPPr, 2);
            break;
        case 193:
            dJydy[0] = (-1.0*myFFr + 1.0*yFFr)/std::pow(sigma_yFFr, 2);
            break;
        case 194:
            dJydy[0] = (-1.0*myEE1E2 + 1.0*yEE1E2)/std::pow(sigma_yEE1E2, 2);
            break;
        case 195:
            dJydy[0] = (-1.0*myEE1Ev3 + 1.0*yEE1Ev3)/std::pow(sigma_yEE1Ev3, 2);
            break;
        case 196:
            dJydy[0] = (-1.0*myEE1E1 + 1.0*yEE1E1)/std::pow(sigma_yEE1E1, 2);
            break;
        case 197:
            dJydy[0] = (-1.0*myEE1E3 + 1.0*yEE1E3)/std::pow(sigma_yEE1E3, 2);
            break;
        case 198:
            dJydy[0] = (-1.0*myEE1E4 + 1.0*yEE1E4)/std::pow(sigma_yEE1E4, 2);
            break;
        case 199:
            dJydy[0] = (-1.0*myE2HE3 + 1.0*yE2HE3)/std::pow(sigma_yE2HE3, 2);
            break;
        case 200:
            dJydy[0] = (-1.0*myE1HE3 + 1.0*yE1HE3)/std::pow(sigma_yE1HE3, 2);
            break;
        case 201:
            dJydy[0] = (-1.0*myHE3E3 + 1.0*yHE3E3)/std::pow(sigma_yHE3E3, 2);
            break;
        case 202:
            dJydy[0] = (-1.0*myHE3Ev3 + 1.0*yHE3Ev3)/std::pow(sigma_yHE3Ev3, 2);
            break;
        case 203:
            dJydy[0] = (-1.0*myHE3E4 + 1.0*yHE3E4)/std::pow(sigma_yHE3E4, 2);
            break;
        case 204:
            dJydy[0] = (-1.0*myE2HE4 + 1.0*yE2HE4)/std::pow(sigma_yE2HE4, 2);
            break;
        case 205:
            dJydy[0] = (-1.0*myHE4Ev3 + 1.0*yHE4Ev3)/std::pow(sigma_yHE4Ev3, 2);
            break;
        case 206:
            dJydy[0] = (-1.0*myE1HE4 + 1.0*yE1HE4)/std::pow(sigma_yE1HE4, 2);
            break;
        case 207:
            dJydy[0] = (-1.0*myE3HE4 + 1.0*yE3HE4)/std::pow(sigma_yE3HE4, 2);
            break;
        case 208:
            dJydy[0] = (-1.0*myHE4E4 + 1.0*yHE4E4)/std::pow(sigma_yHE4E4, 2);
            break;
        case 209:
            dJydy[0] = (-1.0*myHGF_Met_Met + 1.0*yHGF_Met_Met)/std::pow(sigma_yHGF_Met_Met, 2);
            break;
        case 210:
            dJydy[0] = (-1.0*myPPrPr + 1.0*yPPrPr)/std::pow(sigma_yPPrPr, 2);
            break;
        case 211:
            dJydy[0] = (-1.0*myFFrFr + 1.0*yFFrFr)/std::pow(sigma_yFFrFr, 2);
            break;
        case 212:
            dJydy[0] = (-1.0*myIIrIr + 1.0*yIIrIr)/std::pow(sigma_yIIrIr, 2);
            break;
        case 213:
            dJydy[0] = (-1.0*myINS_Isr_Isr + 1.0*yINS_Isr_Isr)/std::pow(sigma_yINS_Isr_Isr, 2);
            break;
        case 214:
            dJydy[0] = (-1.0*myEE1EE1 + 1.0*yEE1EE1)/std::pow(sigma_yEE1EE1, 2);
            break;
        case 215:
            dJydy[0] = (-1.0*myEE1HE3 + 1.0*yEE1HE3)/std::pow(sigma_yEE1HE3, 2);
            break;
        case 216:
            dJydy[0] = (-1.0*myEE1HE4 + 1.0*yEE1HE4)/std::pow(sigma_yEE1HE4, 2);
            break;
        case 217:
            dJydy[0] = (-1.0*myHE3HE3 + 1.0*yHE3HE3)/std::pow(sigma_yHE3HE3, 2);
            break;
        case 218:
            dJydy[0] = (-1.0*myHE3HE4 + 1.0*yHE3HE4)/std::pow(sigma_yHE3HE4, 2);
            break;
        case 219:
            dJydy[0] = (-1.0*myHE4HE4 + 1.0*yHE4HE4)/std::pow(sigma_yHE4HE4, 2);
            break;
        case 220:
            dJydy[0] = (-1.0*myHGF_Met_HGF_Met + 1.0*yHGF_Met_HGF_Met)/std::pow(sigma_yHGF_Met_HGF_Met, 2);
            break;
        case 221:
            dJydy[0] = (-1.0*myPPrPPr + 1.0*yPPrPPr)/std::pow(sigma_yPPrPPr, 2);
            break;
        case 222:
            dJydy[0] = (-1.0*myFFrFFr + 1.0*yFFrFFr)/std::pow(sigma_yFFrFFr, 2);
            break;
        case 223:
            dJydy[0] = (-1.0*myIIrIrI + 1.0*yIIrIrI)/std::pow(sigma_yIIrIrI, 2);
            break;
        case 224:
            dJydy[0] = (-1.0*myINS_Isr_Isr_INS + 1.0*yINS_Isr_Isr_INS)/std::pow(sigma_yINS_Isr_Isr_INS, 2);
            break;
        case 225:
            dJydy[0] = (-1.0*myE1_ppERK + 1.0*yE1_ppERK)/std::pow(sigma_yE1_ppERK, 2);
            break;
        case 226:
            dJydy[0] = (-1.0*myE2_ppERK + 1.0*yE2_ppERK)/std::pow(sigma_yE2_ppERK, 2);
            break;
        case 227:
            dJydy[0] = (-1.0*myE4_ppERK + 1.0*yE4_ppERK)/std::pow(sigma_yE4_ppERK, 2);
            break;
        case 228:
            dJydy[0] = (-1.0*mypEE1E2 + 1.0*ypEE1E2)/std::pow(sigma_ypEE1E2, 2);
            break;
        case 229:
            dJydy[0] = (-1.0*mypEE1Ev3 + 1.0*ypEE1Ev3)/std::pow(sigma_ypEE1Ev3, 2);
            break;
        case 230:
            dJydy[0] = (-1.0*mypEE1E1 + 1.0*ypEE1E1)/std::pow(sigma_ypEE1E1, 2);
            break;
        case 231:
            dJydy[0] = (-1.0*mypEE1EE1 + 1.0*ypEE1EE1)/std::pow(sigma_ypEE1EE1, 2);
            break;
        case 232:
            dJydy[0] = (-1.0*mypEE1E3 + 1.0*ypEE1E3)/std::pow(sigma_ypEE1E3, 2);
            break;
        case 233:
            dJydy[0] = (-1.0*mypEE1HE3 + 1.0*ypEE1HE3)/std::pow(sigma_ypEE1HE3, 2);
            break;
        case 234:
            dJydy[0] = (-1.0*mypEE1E4 + 1.0*ypEE1E4)/std::pow(sigma_ypEE1E4, 2);
            break;
        case 235:
            dJydy[0] = (-1.0*mypEE1HE4 + 1.0*ypEE1HE4)/std::pow(sigma_ypEE1HE4, 2);
            break;
        case 236:
            dJydy[0] = (-1.0*mypE2HE3 + 1.0*ypE2HE3)/std::pow(sigma_ypE2HE3, 2);
            break;
        case 237:
            dJydy[0] = (-1.0*mypHE3Ev3 + 1.0*ypHE3Ev3)/std::pow(sigma_ypHE3Ev3, 2);
            break;
        case 238:
            dJydy[0] = (-1.0*mypE1HE3 + 1.0*ypE1HE3)/std::pow(sigma_ypE1HE3, 2);
            break;
        case 239:
            dJydy[0] = (-1.0*mypHE3E4 + 1.0*ypHE3E4)/std::pow(sigma_ypHE3E4, 2);
            break;
        case 240:
            dJydy[0] = (-1.0*mypHE3HE4 + 1.0*ypHE3HE4)/std::pow(sigma_ypHE3HE4, 2);
            break;
        case 241:
            dJydy[0] = (-1.0*mypE2HE4 + 1.0*ypE2HE4)/std::pow(sigma_ypE2HE4, 2);
            break;
        case 242:
            dJydy[0] = (-1.0*mypHE4Ev3 + 1.0*ypHE4Ev3)/std::pow(sigma_ypHE4Ev3, 2);
            break;
        case 243:
            dJydy[0] = (-1.0*mypE1HE4 + 1.0*ypE1HE4)/std::pow(sigma_ypE1HE4, 2);
            break;
        case 244:
            dJydy[0] = (-1.0*mypE3HE4 + 1.0*ypE3HE4)/std::pow(sigma_ypE3HE4, 2);
            break;
        case 245:
            dJydy[0] = (-1.0*mypHE4E4 + 1.0*ypHE4E4)/std::pow(sigma_ypHE4E4, 2);
            break;
        case 246:
            dJydy[0] = (-1.0*mypHE4HE4 + 1.0*ypHE4HE4)/std::pow(sigma_ypHE4HE4, 2);
            break;
        case 247:
            dJydy[0] = (-1.0*mypHGF_Met_Met + 1.0*ypHGF_Met_Met)/std::pow(sigma_ypHGF_Met_Met, 2);
            break;
        case 248:
            dJydy[0] = (-1.0*mypHGF_Met_HGF_Met + 1.0*ypHGF_Met_HGF_Met)/std::pow(sigma_ypHGF_Met_HGF_Met, 2);
            break;
        case 249:
            dJydy[0] = (-1.0*mypPPrPPr + 1.0*ypPPrPPr)/std::pow(sigma_ypPPrPPr, 2);
            break;
        case 250:
            dJydy[0] = (-1.0*mypPPrPr + 1.0*ypPPrPr)/std::pow(sigma_ypPPrPr, 2);
            break;
        case 251:
            dJydy[0] = (-1.0*mypFFrFFr + 1.0*ypFFrFFr)/std::pow(sigma_ypFFrFFr, 2);
            break;
        case 252:
            dJydy[0] = (-1.0*mypFFrFr + 1.0*ypFFrFr)/std::pow(sigma_ypFFrFr, 2);
            break;
        case 253:
            dJydy[0] = (-1.0*mypIIrIr + 1.0*ypIIrIr)/std::pow(sigma_ypIIrIr, 2);
            break;
        case 254:
            dJydy[0] = (-1.0*mypINS_Isr_Isr + 1.0*ypINS_Isr_Isr)/std::pow(sigma_ypINS_Isr_Isr, 2);
            break;
        case 255:
            dJydy[0] = (-1.0*mypIIrIrI + 1.0*ypIIrIrI)/std::pow(sigma_ypIIrIrI, 2);
            break;
        case 256:
            dJydy[0] = (-1.0*mypINS_Isr_Isr_INS + 1.0*ypINS_Isr_Isr_INS)/std::pow(sigma_ypINS_Isr_Isr_INS, 2);
            break;
        case 257:
            dJydy[0] = (-1.0*mypIIrIr_IRS + 1.0*ypIIrIr_IRS)/std::pow(sigma_ypIIrIr_IRS, 2);
            break;
        case 258:
            dJydy[0] = (-1.0*mypINS_Isr_Isr_IRS + 1.0*ypINS_Isr_Isr_IRS)/std::pow(sigma_ypINS_Isr_Isr_IRS, 2);
            break;
        case 259:
            dJydy[0] = (-1.0*mypIIrIrI_IRS + 1.0*ypIIrIrI_IRS)/std::pow(sigma_ypIIrIrI_IRS, 2);
            break;
        case 260:
            dJydy[0] = (-1.0*mypINS_Isr_Isr_INS_IRS + 1.0*ypINS_Isr_Isr_INS_IRS)/std::pow(sigma_ypINS_Isr_Isr_INS_IRS, 2);
            break;
        case 261:
            dJydy[0] = (-1.0*mySp_EE1E2 + 1.0*ySp_EE1E2)/std::pow(sigma_ySp_EE1E2, 2);
            break;
        case 262:
            dJydy[0] = (-1.0*mySp_EE1Ev3 + 1.0*ySp_EE1Ev3)/std::pow(sigma_ySp_EE1Ev3, 2);
            break;
        case 263:
            dJydy[0] = (-1.0*mySp_EE1E1 + 1.0*ySp_EE1E1)/std::pow(sigma_ySp_EE1E1, 2);
            break;
        case 264:
            dJydy[0] = (-1.0*mySp_EE1EE1 + 1.0*ySp_EE1EE1)/std::pow(sigma_ySp_EE1EE1, 2);
            break;
        case 265:
            dJydy[0] = (-1.0*mySp_EE1E3 + 1.0*ySp_EE1E3)/std::pow(sigma_ySp_EE1E3, 2);
            break;
        case 266:
            dJydy[0] = (-1.0*mySp_EE1HE3 + 1.0*ySp_EE1HE3)/std::pow(sigma_ySp_EE1HE3, 2);
            break;
        case 267:
            dJydy[0] = (-1.0*mySp_EE1E4 + 1.0*ySp_EE1E4)/std::pow(sigma_ySp_EE1E4, 2);
            break;
        case 268:
            dJydy[0] = (-1.0*mySp_EE1HE4 + 1.0*ySp_EE1HE4)/std::pow(sigma_ySp_EE1HE4, 2);
            break;
        case 269:
            dJydy[0] = (-1.0*mySp_E2HE3 + 1.0*ySp_E2HE3)/std::pow(sigma_ySp_E2HE3, 2);
            break;
        case 270:
            dJydy[0] = (-1.0*mySp_HE3Ev3 + 1.0*ySp_HE3Ev3)/std::pow(sigma_ySp_HE3Ev3, 2);
            break;
        case 271:
            dJydy[0] = (-1.0*mySp_E1HE3 + 1.0*ySp_E1HE3)/std::pow(sigma_ySp_E1HE3, 2);
            break;
        case 272:
            dJydy[0] = (-1.0*mySp_HE3E4 + 1.0*ySp_HE3E4)/std::pow(sigma_ySp_HE3E4, 2);
            break;
        case 273:
            dJydy[0] = (-1.0*mySp_HE3HE4 + 1.0*ySp_HE3HE4)/std::pow(sigma_ySp_HE3HE4, 2);
            break;
        case 274:
            dJydy[0] = (-1.0*mySp_E2HE4 + 1.0*ySp_E2HE4)/std::pow(sigma_ySp_E2HE4, 2);
            break;
        case 275:
            dJydy[0] = (-1.0*mySp_HE4Ev3 + 1.0*ySp_HE4Ev3)/std::pow(sigma_ySp_HE4Ev3, 2);
            break;
        case 276:
            dJydy[0] = (-1.0*mySp_E1HE4 + 1.0*ySp_E1HE4)/std::pow(sigma_ySp_E1HE4, 2);
            break;
        case 277:
            dJydy[0] = (-1.0*mySp_E3HE4 + 1.0*ySp_E3HE4)/std::pow(sigma_ySp_E3HE4, 2);
            break;
        case 278:
            dJydy[0] = (-1.0*mySp_HE4E4 + 1.0*ySp_HE4E4)/std::pow(sigma_ySp_HE4E4, 2);
            break;
        case 279:
            dJydy[0] = (-1.0*mySp_HE4HE4 + 1.0*ySp_HE4HE4)/std::pow(sigma_ySp_HE4HE4, 2);
            break;
        case 280:
            dJydy[0] = (-1.0*mySp_HGF_Met_Met + 1.0*ySp_HGF_Met_Met)/std::pow(sigma_ySp_HGF_Met_Met, 2);
            break;
        case 281:
            dJydy[0] = (-1.0*mySp_HGF_Met_HGF_Met + 1.0*ySp_HGF_Met_HGF_Met)/std::pow(sigma_ySp_HGF_Met_HGF_Met, 2);
            break;
        case 282:
            dJydy[0] = (-1.0*mySp_PPrPPr + 1.0*ySp_PPrPPr)/std::pow(sigma_ySp_PPrPPr, 2);
            break;
        case 283:
            dJydy[0] = (-1.0*mySp_PPrPr + 1.0*ySp_PPrPr)/std::pow(sigma_ySp_PPrPr, 2);
            break;
        case 284:
            dJydy[0] = (-1.0*mySp_FFrFFr + 1.0*ySp_FFrFFr)/std::pow(sigma_ySp_FFrFFr, 2);
            break;
        case 285:
            dJydy[0] = (-1.0*mySp_FFrFr + 1.0*ySp_FFrFr)/std::pow(sigma_ySp_FFrFr, 2);
            break;
        case 286:
            dJydy[0] = (-1.0*mySp_IIrIr + 1.0*ySp_IIrIr)/std::pow(sigma_ySp_IIrIr, 2);
            break;
        case 287:
            dJydy[0] = (-1.0*mySp_INS_Isr_Isr + 1.0*ySp_INS_Isr_Isr)/std::pow(sigma_ySp_INS_Isr_Isr, 2);
            break;
        case 288:
            dJydy[0] = (-1.0*mySp_IIrIrI + 1.0*ySp_IIrIrI)/std::pow(sigma_ySp_IIrIrI, 2);
            break;
        case 289:
            dJydy[0] = (-1.0*mySp_INS_Isr_Isr_INS + 1.0*ySp_INS_Isr_Isr_INS)/std::pow(sigma_ySp_INS_Isr_Isr_INS, 2);
            break;
        case 290:
            dJydy[0] = (-1.0*myEE1E2int + 1.0*yEE1E2int)/std::pow(sigma_yEE1E2int, 2);
            break;
        case 291:
            dJydy[0] = (-1.0*myEE1Ev3int + 1.0*yEE1Ev3int)/std::pow(sigma_yEE1Ev3int, 2);
            break;
        case 292:
            dJydy[0] = (-1.0*myEE1E1int + 1.0*yEE1E1int)/std::pow(sigma_yEE1E1int, 2);
            break;
        case 293:
            dJydy[0] = (-1.0*myEE1EE1int + 1.0*yEE1EE1int)/std::pow(sigma_yEE1EE1int, 2);
            break;
        case 294:
            dJydy[0] = (-1.0*myEE1E3int + 1.0*yEE1E3int)/std::pow(sigma_yEE1E3int, 2);
            break;
        case 295:
            dJydy[0] = (-1.0*myEE1HE3int + 1.0*yEE1HE3int)/std::pow(sigma_yEE1HE3int, 2);
            break;
        case 296:
            dJydy[0] = (-1.0*myEE1E4int + 1.0*yEE1E4int)/std::pow(sigma_yEE1E4int, 2);
            break;
        case 297:
            dJydy[0] = (-1.0*myEE1HE4int + 1.0*yEE1HE4int)/std::pow(sigma_yEE1HE4int, 2);
            break;
        case 298:
            dJydy[0] = (-1.0*myE2HE3int + 1.0*yE2HE3int)/std::pow(sigma_yE2HE3int, 2);
            break;
        case 299:
            dJydy[0] = (-1.0*myHE3Ev3int + 1.0*yHE3Ev3int)/std::pow(sigma_yHE3Ev3int, 2);
            break;
        case 300:
            dJydy[0] = (-1.0*myE1HE3int + 1.0*yE1HE3int)/std::pow(sigma_yE1HE3int, 2);
            break;
        case 301:
            dJydy[0] = (-1.0*myHE3E4int + 1.0*yHE3E4int)/std::pow(sigma_yHE3E4int, 2);
            break;
        case 302:
            dJydy[0] = (-1.0*myHE3HE4int + 1.0*yHE3HE4int)/std::pow(sigma_yHE3HE4int, 2);
            break;
        case 303:
            dJydy[0] = (-1.0*myE2HE4int + 1.0*yE2HE4int)/std::pow(sigma_yE2HE4int, 2);
            break;
        case 304:
            dJydy[0] = (-1.0*myHE4Ev3int + 1.0*yHE4Ev3int)/std::pow(sigma_yHE4Ev3int, 2);
            break;
        case 305:
            dJydy[0] = (-1.0*myE1HE4int + 1.0*yE1HE4int)/std::pow(sigma_yE1HE4int, 2);
            break;
        case 306:
            dJydy[0] = (-1.0*myE3HE4int + 1.0*yE3HE4int)/std::pow(sigma_yE3HE4int, 2);
            break;
        case 307:
            dJydy[0] = (-1.0*myHE4E4int + 1.0*yHE4E4int)/std::pow(sigma_yHE4E4int, 2);
            break;
        case 308:
            dJydy[0] = (-1.0*myHE4HE4int + 1.0*yHE4HE4int)/std::pow(sigma_yHE4HE4int, 2);
            break;
        case 309:
            dJydy[0] = (-1.0*myHGF_Met_Metint + 1.0*yHGF_Met_Metint)/std::pow(sigma_yHGF_Met_Metint, 2);
            break;
        case 310:
            dJydy[0] = (-1.0*myHGF_Met_HGF_Metint + 1.0*yHGF_Met_HGF_Metint)/std::pow(sigma_yHGF_Met_HGF_Metint, 2);
            break;
        case 311:
            dJydy[0] = (-1.0*myPPrPPrint + 1.0*yPPrPPrint)/std::pow(sigma_yPPrPPrint, 2);
            break;
        case 312:
            dJydy[0] = (-1.0*myPPrPrint + 1.0*yPPrPrint)/std::pow(sigma_yPPrPrint, 2);
            break;
        case 313:
            dJydy[0] = (-1.0*myFFrFFrint + 1.0*yFFrFFrint)/std::pow(sigma_yFFrFFrint, 2);
            break;
        case 314:
            dJydy[0] = (-1.0*myFFrFrint + 1.0*yFFrFrint)/std::pow(sigma_yFFrFrint, 2);
            break;
        case 315:
            dJydy[0] = (-1.0*myIIrIr_int + 1.0*yIIrIr_int)/std::pow(sigma_yIIrIr_int, 2);
            break;
        case 316:
            dJydy[0] = (-1.0*myINS_Isr_Isr_int + 1.0*yINS_Isr_Isr_int)/std::pow(sigma_yINS_Isr_Isr_int, 2);
            break;
        case 317:
            dJydy[0] = (-1.0*myIIrIrI_int + 1.0*yIIrIrI_int)/std::pow(sigma_yIIrIrI_int, 2);
            break;
        case 318:
            dJydy[0] = (-1.0*myINS_Isr_Isr_INS_int + 1.0*yINS_Isr_Isr_INS_int)/std::pow(sigma_yINS_Isr_Isr_INS_int, 2);
            break;
        case 319:
            dJydy[0] = (-1.0*mypEE1E2int + 1.0*ypEE1E2int)/std::pow(sigma_ypEE1E2int, 2);
            break;
        case 320:
            dJydy[0] = (-1.0*mypEE1Ev3int + 1.0*ypEE1Ev3int)/std::pow(sigma_ypEE1Ev3int, 2);
            break;
        case 321:
            dJydy[0] = (-1.0*mypEE1E1int + 1.0*ypEE1E1int)/std::pow(sigma_ypEE1E1int, 2);
            break;
        case 322:
            dJydy[0] = (-1.0*mypEE1EE1int + 1.0*ypEE1EE1int)/std::pow(sigma_ypEE1EE1int, 2);
            break;
        case 323:
            dJydy[0] = (-1.0*mypEE1E3int + 1.0*ypEE1E3int)/std::pow(sigma_ypEE1E3int, 2);
            break;
        case 324:
            dJydy[0] = (-1.0*mypEE1HE3int + 1.0*ypEE1HE3int)/std::pow(sigma_ypEE1HE3int, 2);
            break;
        case 325:
            dJydy[0] = (-1.0*mypEE1E4int + 1.0*ypEE1E4int)/std::pow(sigma_ypEE1E4int, 2);
            break;
        case 326:
            dJydy[0] = (-1.0*mypEE1HE4int + 1.0*ypEE1HE4int)/std::pow(sigma_ypEE1HE4int, 2);
            break;
        case 327:
            dJydy[0] = (-1.0*mypE2HE3int + 1.0*ypE2HE3int)/std::pow(sigma_ypE2HE3int, 2);
            break;
        case 328:
            dJydy[0] = (-1.0*mypHE3Ev3int + 1.0*ypHE3Ev3int)/std::pow(sigma_ypHE3Ev3int, 2);
            break;
        case 329:
            dJydy[0] = (-1.0*mypE1HE3int + 1.0*ypE1HE3int)/std::pow(sigma_ypE1HE3int, 2);
            break;
        case 330:
            dJydy[0] = (-1.0*mypHE3E4int + 1.0*ypHE3E4int)/std::pow(sigma_ypHE3E4int, 2);
            break;
        case 331:
            dJydy[0] = (-1.0*mypHE3HE4int + 1.0*ypHE3HE4int)/std::pow(sigma_ypHE3HE4int, 2);
            break;
        case 332:
            dJydy[0] = (-1.0*mypE2HE4int + 1.0*ypE2HE4int)/std::pow(sigma_ypE2HE4int, 2);
            break;
        case 333:
            dJydy[0] = (-1.0*mypHE4Ev3int + 1.0*ypHE4Ev3int)/std::pow(sigma_ypHE4Ev3int, 2);
            break;
        case 334:
            dJydy[0] = (-1.0*mypE1HE4int + 1.0*ypE1HE4int)/std::pow(sigma_ypE1HE4int, 2);
            break;
        case 335:
            dJydy[0] = (-1.0*mypE3HE4int + 1.0*ypE3HE4int)/std::pow(sigma_ypE3HE4int, 2);
            break;
        case 336:
            dJydy[0] = (-1.0*mypHE4E4int + 1.0*ypHE4E4int)/std::pow(sigma_ypHE4E4int, 2);
            break;
        case 337:
            dJydy[0] = (-1.0*mypHE4HE4int + 1.0*ypHE4HE4int)/std::pow(sigma_ypHE4HE4int, 2);
            break;
        case 338:
            dJydy[0] = (-1.0*mypHGF_Met_Metint + 1.0*ypHGF_Met_Metint)/std::pow(sigma_ypHGF_Met_Metint, 2);
            break;
        case 339:
            dJydy[0] = (-1.0*mypHGF_Met_HGF_Metint + 1.0*ypHGF_Met_HGF_Metint)/std::pow(sigma_ypHGF_Met_HGF_Metint, 2);
            break;
        case 340:
            dJydy[0] = (-1.0*mypPPrPPrint + 1.0*ypPPrPPrint)/std::pow(sigma_ypPPrPPrint, 2);
            break;
        case 341:
            dJydy[0] = (-1.0*mypPPrPrint + 1.0*ypPPrPrint)/std::pow(sigma_ypPPrPrint, 2);
            break;
        case 342:
            dJydy[0] = (-1.0*mypFFrFFrint + 1.0*ypFFrFFrint)/std::pow(sigma_ypFFrFFrint, 2);
            break;
        case 343:
            dJydy[0] = (-1.0*mypFFrFrint + 1.0*ypFFrFrint)/std::pow(sigma_ypFFrFrint, 2);
            break;
        case 344:
            dJydy[0] = (-1.0*mypIIrIr_int + 1.0*ypIIrIr_int)/std::pow(sigma_ypIIrIr_int, 2);
            break;
        case 345:
            dJydy[0] = (-1.0*mypINS_Isr_Isr_int + 1.0*ypINS_Isr_Isr_int)/std::pow(sigma_ypINS_Isr_Isr_int, 2);
            break;
        case 346:
            dJydy[0] = (-1.0*mypIIrIrI_int + 1.0*ypIIrIrI_int)/std::pow(sigma_ypIIrIrI_int, 2);
            break;
        case 347:
            dJydy[0] = (-1.0*mypINS_Isr_Isr_INS_int + 1.0*ypINS_Isr_Isr_INS_int)/std::pow(sigma_ypINS_Isr_Isr_INS_int, 2);
            break;
        case 348:
            dJydy[0] = (-1.0*mypIIrIr_int_IRS + 1.0*ypIIrIr_int_IRS)/std::pow(sigma_ypIIrIr_int_IRS, 2);
            break;
        case 349:
            dJydy[0] = (-1.0*mypINS_Isr_Isr_int_IRS + 1.0*ypINS_Isr_Isr_int_IRS)/std::pow(sigma_ypINS_Isr_Isr_int_IRS, 2);
            break;
        case 350:
            dJydy[0] = (-1.0*mypIIrIrI_int_IRS + 1.0*ypIIrIrI_int_IRS)/std::pow(sigma_ypIIrIrI_int_IRS, 2);
            break;
        case 351:
            dJydy[0] = (-1.0*mypINS_Isr_Isr_INS_int_IRS + 1.0*ypINS_Isr_Isr_INS_int_IRS)/std::pow(sigma_ypINS_Isr_Isr_INS_int_IRS, 2);
            break;
        case 352:
            dJydy[0] = (-1.0*mypEE1E2_G2_SOS + 1.0*ypEE1E2_G2_SOS)/std::pow(sigma_ypEE1E2_G2_SOS, 2);
            break;
        case 353:
            dJydy[0] = (-1.0*mypEE1Ev3_G2_SOS + 1.0*ypEE1Ev3_G2_SOS)/std::pow(sigma_ypEE1Ev3_G2_SOS, 2);
            break;
        case 354:
            dJydy[0] = (-1.0*mypEE1E1_G2_SOS + 1.0*ypEE1E1_G2_SOS)/std::pow(sigma_ypEE1E1_G2_SOS, 2);
            break;
        case 355:
            dJydy[0] = (-1.0*mypEE1EE1_G2_SOS + 1.0*ypEE1EE1_G2_SOS)/std::pow(sigma_ypEE1EE1_G2_SOS, 2);
            break;
        case 356:
            dJydy[0] = (-1.0*mypEE1E3_G2_SOS + 1.0*ypEE1E3_G2_SOS)/std::pow(sigma_ypEE1E3_G2_SOS, 2);
            break;
        case 357:
            dJydy[0] = (-1.0*mypEE1HE3_G2_SOS + 1.0*ypEE1HE3_G2_SOS)/std::pow(sigma_ypEE1HE3_G2_SOS, 2);
            break;
        case 358:
            dJydy[0] = (-1.0*mypEE1E4_G2_SOS + 1.0*ypEE1E4_G2_SOS)/std::pow(sigma_ypEE1E4_G2_SOS, 2);
            break;
        case 359:
            dJydy[0] = (-1.0*mypEE1HE4_G2_SOS + 1.0*ypEE1HE4_G2_SOS)/std::pow(sigma_ypEE1HE4_G2_SOS, 2);
            break;
        case 360:
            dJydy[0] = (-1.0*mypE2HE3_G2_SOS + 1.0*ypE2HE3_G2_SOS)/std::pow(sigma_ypE2HE3_G2_SOS, 2);
            break;
        case 361:
            dJydy[0] = (-1.0*mypHE3Ev3_G2_SOS + 1.0*ypHE3Ev3_G2_SOS)/std::pow(sigma_ypHE3Ev3_G2_SOS, 2);
            break;
        case 362:
            dJydy[0] = (-1.0*mypE1HE3_G2_SOS + 1.0*ypE1HE3_G2_SOS)/std::pow(sigma_ypE1HE3_G2_SOS, 2);
            break;
        case 363:
            dJydy[0] = (-1.0*mypHE3E4_G2_SOS + 1.0*ypHE3E4_G2_SOS)/std::pow(sigma_ypHE3E4_G2_SOS, 2);
            break;
        case 364:
            dJydy[0] = (-1.0*mypHE3HE4_G2_SOS + 1.0*ypHE3HE4_G2_SOS)/std::pow(sigma_ypHE3HE4_G2_SOS, 2);
            break;
        case 365:
            dJydy[0] = (-1.0*mypE2HE4_G2_SOS + 1.0*ypE2HE4_G2_SOS)/std::pow(sigma_ypE2HE4_G2_SOS, 2);
            break;
        case 366:
            dJydy[0] = (-1.0*mypHE4Ev3_G2_SOS + 1.0*ypHE4Ev3_G2_SOS)/std::pow(sigma_ypHE4Ev3_G2_SOS, 2);
            break;
        case 367:
            dJydy[0] = (-1.0*mypE1HE4_G2_SOS + 1.0*ypE1HE4_G2_SOS)/std::pow(sigma_ypE1HE4_G2_SOS, 2);
            break;
        case 368:
            dJydy[0] = (-1.0*mypE3HE4_G2_SOS + 1.0*ypE3HE4_G2_SOS)/std::pow(sigma_ypE3HE4_G2_SOS, 2);
            break;
        case 369:
            dJydy[0] = (-1.0*mypHE4E4_G2_SOS + 1.0*ypHE4E4_G2_SOS)/std::pow(sigma_ypHE4E4_G2_SOS, 2);
            break;
        case 370:
            dJydy[0] = (-1.0*mypHE4HE4_G2_SOS + 1.0*ypHE4HE4_G2_SOS)/std::pow(sigma_ypHE4HE4_G2_SOS, 2);
            break;
        case 371:
            dJydy[0] = (-1.0*mypHGF_Met_Met_G2_SOS + 1.0*ypHGF_Met_Met_G2_SOS)/std::pow(sigma_ypHGF_Met_Met_G2_SOS, 2);
            break;
        case 372:
            dJydy[0] = (-1.0*mypHGF_Met_HGF_Met_G2_SOS + 1.0*ypHGF_Met_HGF_Met_G2_SOS)/std::pow(sigma_ypHGF_Met_HGF_Met_G2_SOS, 2);
            break;
        case 373:
            dJydy[0] = (-1.0*mypPPrPPr_G2_SOS + 1.0*ypPPrPPr_G2_SOS)/std::pow(sigma_ypPPrPPr_G2_SOS, 2);
            break;
        case 374:
            dJydy[0] = (-1.0*mypPPrPr_G2_SOS + 1.0*ypPPrPr_G2_SOS)/std::pow(sigma_ypPPrPr_G2_SOS, 2);
            break;
        case 375:
            dJydy[0] = (-1.0*mypFFrFFr_G2_SOS + 1.0*ypFFrFFr_G2_SOS)/std::pow(sigma_ypFFrFFr_G2_SOS, 2);
            break;
        case 376:
            dJydy[0] = (-1.0*mypFFrFr_G2_SOS + 1.0*ypFFrFr_G2_SOS)/std::pow(sigma_ypFFrFr_G2_SOS, 2);
            break;
        case 377:
            dJydy[0] = (-1.0*mypIIrIr_IRS_G2_SOS + 1.0*ypIIrIr_IRS_G2_SOS)/std::pow(sigma_ypIIrIr_IRS_G2_SOS, 2);
            break;
        case 378:
            dJydy[0] = (-1.0*mypINS_Isr_Isr_IRS_G2_SOS + 1.0*ypINS_Isr_Isr_IRS_G2_SOS)/std::pow(sigma_ypINS_Isr_Isr_IRS_G2_SOS, 2);
            break;
        case 379:
            dJydy[0] = (-1.0*mypIIrIrI_IRS_G2_SOS + 1.0*ypIIrIrI_IRS_G2_SOS)/std::pow(sigma_ypIIrIrI_IRS_G2_SOS, 2);
            break;
        case 380:
            dJydy[0] = (-1.0*mypINS_Isr_Isr_INS_IRS_G2_SOS + 1.0*ypINS_Isr_Isr_INS_IRS_G2_SOS)/std::pow(sigma_ypINS_Isr_Isr_INS_IRS_G2_SOS, 2);
            break;
        case 381:
            dJydy[0] = (-1.0*mypEE1E2int_G2_SOS + 1.0*ypEE1E2int_G2_SOS)/std::pow(sigma_ypEE1E2int_G2_SOS, 2);
            break;
        case 382:
            dJydy[0] = (-1.0*mypEE1Ev3int_G2_SOS + 1.0*ypEE1Ev3int_G2_SOS)/std::pow(sigma_ypEE1Ev3int_G2_SOS, 2);
            break;
        case 383:
            dJydy[0] = (-1.0*mypEE1E1int_G2_SOS + 1.0*ypEE1E1int_G2_SOS)/std::pow(sigma_ypEE1E1int_G2_SOS, 2);
            break;
        case 384:
            dJydy[0] = (-1.0*mypEE1EE1int_G2_SOS + 1.0*ypEE1EE1int_G2_SOS)/std::pow(sigma_ypEE1EE1int_G2_SOS, 2);
            break;
        case 385:
            dJydy[0] = (-1.0*mypEE1E3int_G2_SOS + 1.0*ypEE1E3int_G2_SOS)/std::pow(sigma_ypEE1E3int_G2_SOS, 2);
            break;
        case 386:
            dJydy[0] = (-1.0*mypEE1HE3int_G2_SOS + 1.0*ypEE1HE3int_G2_SOS)/std::pow(sigma_ypEE1HE3int_G2_SOS, 2);
            break;
        case 387:
            dJydy[0] = (-1.0*mypEE1E4int_G2_SOS + 1.0*ypEE1E4int_G2_SOS)/std::pow(sigma_ypEE1E4int_G2_SOS, 2);
            break;
        case 388:
            dJydy[0] = (-1.0*mypEE1HE4int_G2_SOS + 1.0*ypEE1HE4int_G2_SOS)/std::pow(sigma_ypEE1HE4int_G2_SOS, 2);
            break;
        case 389:
            dJydy[0] = (-1.0*mypE2HE3int_G2_SOS + 1.0*ypE2HE3int_G2_SOS)/std::pow(sigma_ypE2HE3int_G2_SOS, 2);
            break;
        case 390:
            dJydy[0] = (-1.0*mypHE3Ev3int_G2_SOS + 1.0*ypHE3Ev3int_G2_SOS)/std::pow(sigma_ypHE3Ev3int_G2_SOS, 2);
            break;
        case 391:
            dJydy[0] = (-1.0*mypE1HE3int_G2_SOS + 1.0*ypE1HE3int_G2_SOS)/std::pow(sigma_ypE1HE3int_G2_SOS, 2);
            break;
        case 392:
            dJydy[0] = (-1.0*mypHE3E4int_G2_SOS + 1.0*ypHE3E4int_G2_SOS)/std::pow(sigma_ypHE3E4int_G2_SOS, 2);
            break;
        case 393:
            dJydy[0] = (-1.0*mypHE3HE4int_G2_SOS + 1.0*ypHE3HE4int_G2_SOS)/std::pow(sigma_ypHE3HE4int_G2_SOS, 2);
            break;
        case 394:
            dJydy[0] = (-1.0*mypE2HE4int_G2_SOS + 1.0*ypE2HE4int_G2_SOS)/std::pow(sigma_ypE2HE4int_G2_SOS, 2);
            break;
        case 395:
            dJydy[0] = (-1.0*mypHE4Ev3int_G2_SOS + 1.0*ypHE4Ev3int_G2_SOS)/std::pow(sigma_ypHE4Ev3int_G2_SOS, 2);
            break;
        case 396:
            dJydy[0] = (-1.0*mypE1HE4int_G2_SOS + 1.0*ypE1HE4int_G2_SOS)/std::pow(sigma_ypE1HE4int_G2_SOS, 2);
            break;
        case 397:
            dJydy[0] = (-1.0*mypE3HE4int_G2_SOS + 1.0*ypE3HE4int_G2_SOS)/std::pow(sigma_ypE3HE4int_G2_SOS, 2);
            break;
        case 398:
            dJydy[0] = (-1.0*mypHE4E4int_G2_SOS + 1.0*ypHE4E4int_G2_SOS)/std::pow(sigma_ypHE4E4int_G2_SOS, 2);
            break;
        case 399:
            dJydy[0] = (-1.0*mypHE4HE4int_G2_SOS + 1.0*ypHE4HE4int_G2_SOS)/std::pow(sigma_ypHE4HE4int_G2_SOS, 2);
            break;
        case 400:
            dJydy[0] = (-1.0*mypHGF_Met_Metint_G2_SOS + 1.0*ypHGF_Met_Metint_G2_SOS)/std::pow(sigma_ypHGF_Met_Metint_G2_SOS, 2);
            break;
        case 401:
            dJydy[0] = (-1.0*mypHGF_Met_HGF_Metint_G2_SOS + 1.0*ypHGF_Met_HGF_Metint_G2_SOS)/std::pow(sigma_ypHGF_Met_HGF_Metint_G2_SOS, 2);
            break;
        case 402:
            dJydy[0] = (-1.0*mypPPrPPrint_G2_SOS + 1.0*ypPPrPPrint_G2_SOS)/std::pow(sigma_ypPPrPPrint_G2_SOS, 2);
            break;
        case 403:
            dJydy[0] = (-1.0*mypPPrPrint_G2_SOS + 1.0*ypPPrPrint_G2_SOS)/std::pow(sigma_ypPPrPrint_G2_SOS, 2);
            break;
        case 404:
            dJydy[0] = (-1.0*mypFFrFFrint_G2_SOS + 1.0*ypFFrFFrint_G2_SOS)/std::pow(sigma_ypFFrFFrint_G2_SOS, 2);
            break;
        case 405:
            dJydy[0] = (-1.0*mypFFrFrint_G2_SOS + 1.0*ypFFrFrint_G2_SOS)/std::pow(sigma_ypFFrFrint_G2_SOS, 2);
            break;
        case 406:
            dJydy[0] = (-1.0*mypIIrIr_int_IRS_G2_SOS + 1.0*ypIIrIr_int_IRS_G2_SOS)/std::pow(sigma_ypIIrIr_int_IRS_G2_SOS, 2);
            break;
        case 407:
            dJydy[0] = (-1.0*mypINS_Isr_Isr_int_IRS_G2_SOS + 1.0*ypINS_Isr_Isr_int_IRS_G2_SOS)/std::pow(sigma_ypINS_Isr_Isr_int_IRS_G2_SOS, 2);
            break;
        case 408:
            dJydy[0] = (-1.0*mypIIrIrI_int_IRS_G2_SOS + 1.0*ypIIrIrI_int_IRS_G2_SOS)/std::pow(sigma_ypIIrIrI_int_IRS_G2_SOS, 2);
            break;
        case 409:
            dJydy[0] = (-1.0*mypINS_Isr_Isr_INS_int_IRS_G2_SOS + 1.0*ypINS_Isr_Isr_INS_int_IRS_G2_SOS)/std::pow(sigma_ypINS_Isr_Isr_INS_int_IRS_G2_SOS, 2);
            break;
        case 410:
            dJydy[0] = (-1.0*mypEE1E2_PLCg + 1.0*ypEE1E2_PLCg)/std::pow(sigma_ypEE1E2_PLCg, 2);
            break;
        case 411:
            dJydy[0] = (-1.0*mypEE1Ev3_PLCg + 1.0*ypEE1Ev3_PLCg)/std::pow(sigma_ypEE1Ev3_PLCg, 2);
            break;
        case 412:
            dJydy[0] = (-1.0*mypEE1E1_PLCg + 1.0*ypEE1E1_PLCg)/std::pow(sigma_ypEE1E1_PLCg, 2);
            break;
        case 413:
            dJydy[0] = (-1.0*mypEE1EE1_PLCg + 1.0*ypEE1EE1_PLCg)/std::pow(sigma_ypEE1EE1_PLCg, 2);
            break;
        case 414:
            dJydy[0] = (-1.0*mypEE1E3_PLCg + 1.0*ypEE1E3_PLCg)/std::pow(sigma_ypEE1E3_PLCg, 2);
            break;
        case 415:
            dJydy[0] = (-1.0*mypEE1HE3_PLCg + 1.0*ypEE1HE3_PLCg)/std::pow(sigma_ypEE1HE3_PLCg, 2);
            break;
        case 416:
            dJydy[0] = (-1.0*mypEE1E4_PLCg + 1.0*ypEE1E4_PLCg)/std::pow(sigma_ypEE1E4_PLCg, 2);
            break;
        case 417:
            dJydy[0] = (-1.0*mypEE1HE4_PLCg + 1.0*ypEE1HE4_PLCg)/std::pow(sigma_ypEE1HE4_PLCg, 2);
            break;
        case 418:
            dJydy[0] = (-1.0*mypE2HE3_PLCg + 1.0*ypE2HE3_PLCg)/std::pow(sigma_ypE2HE3_PLCg, 2);
            break;
        case 419:
            dJydy[0] = (-1.0*mypHE3Ev3_PLCg + 1.0*ypHE3Ev3_PLCg)/std::pow(sigma_ypHE3Ev3_PLCg, 2);
            break;
        case 420:
            dJydy[0] = (-1.0*mypE1HE3_PLCg + 1.0*ypE1HE3_PLCg)/std::pow(sigma_ypE1HE3_PLCg, 2);
            break;
        case 421:
            dJydy[0] = (-1.0*mypHE3E4_PLCg + 1.0*ypHE3E4_PLCg)/std::pow(sigma_ypHE3E4_PLCg, 2);
            break;
        case 422:
            dJydy[0] = (-1.0*mypHE3HE4_PLCg + 1.0*ypHE3HE4_PLCg)/std::pow(sigma_ypHE3HE4_PLCg, 2);
            break;
        case 423:
            dJydy[0] = (-1.0*mypE2HE4_PLCg + 1.0*ypE2HE4_PLCg)/std::pow(sigma_ypE2HE4_PLCg, 2);
            break;
        case 424:
            dJydy[0] = (-1.0*mypHE4Ev3_PLCg + 1.0*ypHE4Ev3_PLCg)/std::pow(sigma_ypHE4Ev3_PLCg, 2);
            break;
        case 425:
            dJydy[0] = (-1.0*mypE1HE4_PLCg + 1.0*ypE1HE4_PLCg)/std::pow(sigma_ypE1HE4_PLCg, 2);
            break;
        case 426:
            dJydy[0] = (-1.0*mypE3HE4_PLCg + 1.0*ypE3HE4_PLCg)/std::pow(sigma_ypE3HE4_PLCg, 2);
            break;
        case 427:
            dJydy[0] = (-1.0*mypHE4E4_PLCg + 1.0*ypHE4E4_PLCg)/std::pow(sigma_ypHE4E4_PLCg, 2);
            break;
        case 428:
            dJydy[0] = (-1.0*mypHE4HE4_PLCg + 1.0*ypHE4HE4_PLCg)/std::pow(sigma_ypHE4HE4_PLCg, 2);
            break;
        case 429:
            dJydy[0] = (-1.0*mypHGF_Met_Met_PLCg + 1.0*ypHGF_Met_Met_PLCg)/std::pow(sigma_ypHGF_Met_Met_PLCg, 2);
            break;
        case 430:
            dJydy[0] = (-1.0*mypHGF_Met_HGF_Met_PLCg + 1.0*ypHGF_Met_HGF_Met_PLCg)/std::pow(sigma_ypHGF_Met_HGF_Met_PLCg, 2);
            break;
        case 431:
            dJydy[0] = (-1.0*mypPPrPPr_PLCg + 1.0*ypPPrPPr_PLCg)/std::pow(sigma_ypPPrPPr_PLCg, 2);
            break;
        case 432:
            dJydy[0] = (-1.0*mypPPrPr_PLCg + 1.0*ypPPrPr_PLCg)/std::pow(sigma_ypPPrPr_PLCg, 2);
            break;
        case 433:
            dJydy[0] = (-1.0*mypFFrFFr_PLCg + 1.0*ypFFrFFr_PLCg)/std::pow(sigma_ypFFrFFr_PLCg, 2);
            break;
        case 434:
            dJydy[0] = (-1.0*mypFFrFr_PLCg + 1.0*ypFFrFr_PLCg)/std::pow(sigma_ypFFrFr_PLCg, 2);
            break;
        case 435:
            dJydy[0] = (-1.0*mypIIrIr_IRS_PLCg + 1.0*ypIIrIr_IRS_PLCg)/std::pow(sigma_ypIIrIr_IRS_PLCg, 2);
            break;
        case 436:
            dJydy[0] = (-1.0*mypINS_Isr_Isr_IRS_PLCg + 1.0*ypINS_Isr_Isr_IRS_PLCg)/std::pow(sigma_ypINS_Isr_Isr_IRS_PLCg, 2);
            break;
        case 437:
            dJydy[0] = (-1.0*mypIIrIrI_IRS_PLCg + 1.0*ypIIrIrI_IRS_PLCg)/std::pow(sigma_ypIIrIrI_IRS_PLCg, 2);
            break;
        case 438:
            dJydy[0] = (-1.0*mypINS_Isr_Isr_INS_IRS_PLCg + 1.0*ypINS_Isr_Isr_INS_IRS_PLCg)/std::pow(sigma_ypINS_Isr_Isr_INS_IRS_PLCg, 2);
            break;
        case 439:
            dJydy[0] = (-1.0*mypEE1E2_PI3K1 + 1.0*ypEE1E2_PI3K1)/std::pow(sigma_ypEE1E2_PI3K1, 2);
            break;
        case 440:
            dJydy[0] = (-1.0*mypEE1Ev3_PI3K1 + 1.0*ypEE1Ev3_PI3K1)/std::pow(sigma_ypEE1Ev3_PI3K1, 2);
            break;
        case 441:
            dJydy[0] = (-1.0*mypEE1E1_PI3K1 + 1.0*ypEE1E1_PI3K1)/std::pow(sigma_ypEE1E1_PI3K1, 2);
            break;
        case 442:
            dJydy[0] = (-1.0*mypEE1EE1_PI3K1 + 1.0*ypEE1EE1_PI3K1)/std::pow(sigma_ypEE1EE1_PI3K1, 2);
            break;
        case 443:
            dJydy[0] = (-1.0*mypEE1E3_PI3K1 + 1.0*ypEE1E3_PI3K1)/std::pow(sigma_ypEE1E3_PI3K1, 2);
            break;
        case 444:
            dJydy[0] = (-1.0*mypEE1HE3_PI3K1 + 1.0*ypEE1HE3_PI3K1)/std::pow(sigma_ypEE1HE3_PI3K1, 2);
            break;
        case 445:
            dJydy[0] = (-1.0*mypEE1E4_PI3K1 + 1.0*ypEE1E4_PI3K1)/std::pow(sigma_ypEE1E4_PI3K1, 2);
            break;
        case 446:
            dJydy[0] = (-1.0*mypEE1HE4_PI3K1 + 1.0*ypEE1HE4_PI3K1)/std::pow(sigma_ypEE1HE4_PI3K1, 2);
            break;
        case 447:
            dJydy[0] = (-1.0*mypE2HE3_PI3K1 + 1.0*ypE2HE3_PI3K1)/std::pow(sigma_ypE2HE3_PI3K1, 2);
            break;
        case 448:
            dJydy[0] = (-1.0*mypHE3Ev3_PI3K1 + 1.0*ypHE3Ev3_PI3K1)/std::pow(sigma_ypHE3Ev3_PI3K1, 2);
            break;
        case 449:
            dJydy[0] = (-1.0*mypE1HE3_PI3K1 + 1.0*ypE1HE3_PI3K1)/std::pow(sigma_ypE1HE3_PI3K1, 2);
            break;
        case 450:
            dJydy[0] = (-1.0*mypHE3E4_PI3K1 + 1.0*ypHE3E4_PI3K1)/std::pow(sigma_ypHE3E4_PI3K1, 2);
            break;
        case 451:
            dJydy[0] = (-1.0*mypHE3HE4_PI3K1 + 1.0*ypHE3HE4_PI3K1)/std::pow(sigma_ypHE3HE4_PI3K1, 2);
            break;
        case 452:
            dJydy[0] = (-1.0*mypE2HE4_PI3K1 + 1.0*ypE2HE4_PI3K1)/std::pow(sigma_ypE2HE4_PI3K1, 2);
            break;
        case 453:
            dJydy[0] = (-1.0*mypHE4Ev3_PI3K1 + 1.0*ypHE4Ev3_PI3K1)/std::pow(sigma_ypHE4Ev3_PI3K1, 2);
            break;
        case 454:
            dJydy[0] = (-1.0*mypE1HE4_PI3K1 + 1.0*ypE1HE4_PI3K1)/std::pow(sigma_ypE1HE4_PI3K1, 2);
            break;
        case 455:
            dJydy[0] = (-1.0*mypE3HE4_PI3K1 + 1.0*ypE3HE4_PI3K1)/std::pow(sigma_ypE3HE4_PI3K1, 2);
            break;
        case 456:
            dJydy[0] = (-1.0*mypHE4E4_PI3K1 + 1.0*ypHE4E4_PI3K1)/std::pow(sigma_ypHE4E4_PI3K1, 2);
            break;
        case 457:
            dJydy[0] = (-1.0*mypHE4HE4_PI3K1 + 1.0*ypHE4HE4_PI3K1)/std::pow(sigma_ypHE4HE4_PI3K1, 2);
            break;
        case 458:
            dJydy[0] = (-1.0*mypHGF_Met_Met_PI3K1 + 1.0*ypHGF_Met_Met_PI3K1)/std::pow(sigma_ypHGF_Met_Met_PI3K1, 2);
            break;
        case 459:
            dJydy[0] = (-1.0*mypHGF_Met_HGF_Met_PI3K1 + 1.0*ypHGF_Met_HGF_Met_PI3K1)/std::pow(sigma_ypHGF_Met_HGF_Met_PI3K1, 2);
            break;
        case 460:
            dJydy[0] = (-1.0*mypPPrPPr_PI3K1 + 1.0*ypPPrPPr_PI3K1)/std::pow(sigma_ypPPrPPr_PI3K1, 2);
            break;
        case 461:
            dJydy[0] = (-1.0*mypPPrPr_PI3K1 + 1.0*ypPPrPr_PI3K1)/std::pow(sigma_ypPPrPr_PI3K1, 2);
            break;
        case 462:
            dJydy[0] = (-1.0*mypFFrFFr_PI3K1 + 1.0*ypFFrFFr_PI3K1)/std::pow(sigma_ypFFrFFr_PI3K1, 2);
            break;
        case 463:
            dJydy[0] = (-1.0*mypFFrFr_PI3K1 + 1.0*ypFFrFr_PI3K1)/std::pow(sigma_ypFFrFr_PI3K1, 2);
            break;
        case 464:
            dJydy[0] = (-1.0*mypIIrIr_IRS_PI3K1 + 1.0*ypIIrIr_IRS_PI3K1)/std::pow(sigma_ypIIrIr_IRS_PI3K1, 2);
            break;
        case 465:
            dJydy[0] = (-1.0*mypINS_Isr_Isr_IRS_PI3K1 + 1.0*ypINS_Isr_Isr_IRS_PI3K1)/std::pow(sigma_ypINS_Isr_Isr_IRS_PI3K1, 2);
            break;
        case 466:
            dJydy[0] = (-1.0*mypIIrIrI_IRS_PI3K1 + 1.0*ypIIrIrI_IRS_PI3K1)/std::pow(sigma_ypIIrIrI_IRS_PI3K1, 2);
            break;
        case 467:
            dJydy[0] = (-1.0*mypINS_Isr_Isr_INS_IRS_PI3K1 + 1.0*ypINS_Isr_Isr_INS_IRS_PI3K1)/std::pow(sigma_ypINS_Isr_Isr_INS_IRS_PI3K1, 2);
            break;
        case 468:
            dJydy[0] = (-1.0*mypEE1E2_PI3K2 + 1.0*ypEE1E2_PI3K2)/std::pow(sigma_ypEE1E2_PI3K2, 2);
            break;
        case 469:
            dJydy[0] = (-1.0*mypEE1Ev3_PI3K2 + 1.0*ypEE1Ev3_PI3K2)/std::pow(sigma_ypEE1Ev3_PI3K2, 2);
            break;
        case 470:
            dJydy[0] = (-1.0*mypEE1E1_PI3K2 + 1.0*ypEE1E1_PI3K2)/std::pow(sigma_ypEE1E1_PI3K2, 2);
            break;
        case 471:
            dJydy[0] = (-1.0*mypEE1EE1_PI3K2 + 1.0*ypEE1EE1_PI3K2)/std::pow(sigma_ypEE1EE1_PI3K2, 2);
            break;
        case 472:
            dJydy[0] = (-1.0*mypEE1E3_PI3K2 + 1.0*ypEE1E3_PI3K2)/std::pow(sigma_ypEE1E3_PI3K2, 2);
            break;
        case 473:
            dJydy[0] = (-1.0*mypEE1HE3_PI3K2 + 1.0*ypEE1HE3_PI3K2)/std::pow(sigma_ypEE1HE3_PI3K2, 2);
            break;
        case 474:
            dJydy[0] = (-1.0*mypEE1E4_PI3K2 + 1.0*ypEE1E4_PI3K2)/std::pow(sigma_ypEE1E4_PI3K2, 2);
            break;
        case 475:
            dJydy[0] = (-1.0*mypEE1HE4_PI3K2 + 1.0*ypEE1HE4_PI3K2)/std::pow(sigma_ypEE1HE4_PI3K2, 2);
            break;
        case 476:
            dJydy[0] = (-1.0*mypE2HE3_PI3K2 + 1.0*ypE2HE3_PI3K2)/std::pow(sigma_ypE2HE3_PI3K2, 2);
            break;
        case 477:
            dJydy[0] = (-1.0*mypHE3Ev3_PI3K2 + 1.0*ypHE3Ev3_PI3K2)/std::pow(sigma_ypHE3Ev3_PI3K2, 2);
            break;
        case 478:
            dJydy[0] = (-1.0*mypE1HE3_PI3K2 + 1.0*ypE1HE3_PI3K2)/std::pow(sigma_ypE1HE3_PI3K2, 2);
            break;
        case 479:
            dJydy[0] = (-1.0*mypHE3E4_PI3K2 + 1.0*ypHE3E4_PI3K2)/std::pow(sigma_ypHE3E4_PI3K2, 2);
            break;
        case 480:
            dJydy[0] = (-1.0*mypHE3HE4_PI3K2 + 1.0*ypHE3HE4_PI3K2)/std::pow(sigma_ypHE3HE4_PI3K2, 2);
            break;
        case 481:
            dJydy[0] = (-1.0*mypE2HE4_PI3K2 + 1.0*ypE2HE4_PI3K2)/std::pow(sigma_ypE2HE4_PI3K2, 2);
            break;
        case 482:
            dJydy[0] = (-1.0*mypHE4Ev3_PI3K2 + 1.0*ypHE4Ev3_PI3K2)/std::pow(sigma_ypHE4Ev3_PI3K2, 2);
            break;
        case 483:
            dJydy[0] = (-1.0*mypE1HE4_PI3K2 + 1.0*ypE1HE4_PI3K2)/std::pow(sigma_ypE1HE4_PI3K2, 2);
            break;
        case 484:
            dJydy[0] = (-1.0*mypE3HE4_PI3K2 + 1.0*ypE3HE4_PI3K2)/std::pow(sigma_ypE3HE4_PI3K2, 2);
            break;
        case 485:
            dJydy[0] = (-1.0*mypHE4E4_PI3K2 + 1.0*ypHE4E4_PI3K2)/std::pow(sigma_ypHE4E4_PI3K2, 2);
            break;
        case 486:
            dJydy[0] = (-1.0*mypHE4HE4_PI3K2 + 1.0*ypHE4HE4_PI3K2)/std::pow(sigma_ypHE4HE4_PI3K2, 2);
            break;
        case 487:
            dJydy[0] = (-1.0*mypHGF_Met_Met_PI3K2 + 1.0*ypHGF_Met_Met_PI3K2)/std::pow(sigma_ypHGF_Met_Met_PI3K2, 2);
            break;
        case 488:
            dJydy[0] = (-1.0*mypHGF_Met_HGF_Met_PI3K2 + 1.0*ypHGF_Met_HGF_Met_PI3K2)/std::pow(sigma_ypHGF_Met_HGF_Met_PI3K2, 2);
            break;
        case 489:
            dJydy[0] = (-1.0*mypPPrPPr_PI3K2 + 1.0*ypPPrPPr_PI3K2)/std::pow(sigma_ypPPrPPr_PI3K2, 2);
            break;
        case 490:
            dJydy[0] = (-1.0*mypPPrPr_PI3K2 + 1.0*ypPPrPr_PI3K2)/std::pow(sigma_ypPPrPr_PI3K2, 2);
            break;
        case 491:
            dJydy[0] = (-1.0*mypFFrFFr_PI3K2 + 1.0*ypFFrFFr_PI3K2)/std::pow(sigma_ypFFrFFr_PI3K2, 2);
            break;
        case 492:
            dJydy[0] = (-1.0*mypFFrFr_PI3K2 + 1.0*ypFFrFr_PI3K2)/std::pow(sigma_ypFFrFr_PI3K2, 2);
            break;
        case 493:
            dJydy[0] = (-1.0*mypIIrIr_IRS_PI3K2 + 1.0*ypIIrIr_IRS_PI3K2)/std::pow(sigma_ypIIrIr_IRS_PI3K2, 2);
            break;
        case 494:
            dJydy[0] = (-1.0*mypINS_Isr_Isr_IRS_PI3K2 + 1.0*ypINS_Isr_Isr_IRS_PI3K2)/std::pow(sigma_ypINS_Isr_Isr_IRS_PI3K2, 2);
            break;
        case 495:
            dJydy[0] = (-1.0*mypIIrIrI_IRS_PI3K2 + 1.0*ypIIrIrI_IRS_PI3K2)/std::pow(sigma_ypIIrIrI_IRS_PI3K2, 2);
            break;
        case 496:
            dJydy[0] = (-1.0*mypINS_Isr_Isr_INS_IRS_PI3K2 + 1.0*ypINS_Isr_Isr_INS_IRS_PI3K2)/std::pow(sigma_ypINS_Isr_Isr_INS_IRS_PI3K2, 2);
            break;
        case 497:
            dJydy[0] = (-1.0*mypEE1E2int_G2_SOS_RasD + 1.0*ypEE1E2int_G2_SOS_RasD)/std::pow(sigma_ypEE1E2int_G2_SOS_RasD, 2);
            break;
        case 498:
            dJydy[0] = (-1.0*mypEE1Ev3int_G2_SOS_RasD + 1.0*ypEE1Ev3int_G2_SOS_RasD)/std::pow(sigma_ypEE1Ev3int_G2_SOS_RasD, 2);
            break;
        case 499:
            dJydy[0] = (-1.0*mypEE1E1int_G2_SOS_RasD + 1.0*ypEE1E1int_G2_SOS_RasD)/std::pow(sigma_ypEE1E1int_G2_SOS_RasD, 2);
            break;
        case 500:
            dJydy[0] = (-1.0*mypEE1EE1int_G2_SOS_RasD + 1.0*ypEE1EE1int_G2_SOS_RasD)/std::pow(sigma_ypEE1EE1int_G2_SOS_RasD, 2);
            break;
        case 501:
            dJydy[0] = (-1.0*mypEE1E3int_G2_SOS_RasD + 1.0*ypEE1E3int_G2_SOS_RasD)/std::pow(sigma_ypEE1E3int_G2_SOS_RasD, 2);
            break;
        case 502:
            dJydy[0] = (-1.0*mypEE1HE3int_G2_SOS_RasD + 1.0*ypEE1HE3int_G2_SOS_RasD)/std::pow(sigma_ypEE1HE3int_G2_SOS_RasD, 2);
            break;
        case 503:
            dJydy[0] = (-1.0*mypEE1E4int_G2_SOS_RasD + 1.0*ypEE1E4int_G2_SOS_RasD)/std::pow(sigma_ypEE1E4int_G2_SOS_RasD, 2);
            break;
        case 504:
            dJydy[0] = (-1.0*mypEE1HE4int_G2_SOS_RasD + 1.0*ypEE1HE4int_G2_SOS_RasD)/std::pow(sigma_ypEE1HE4int_G2_SOS_RasD, 2);
            break;
        case 505:
            dJydy[0] = (-1.0*mypE2HE3int_G2_SOS_RasD + 1.0*ypE2HE3int_G2_SOS_RasD)/std::pow(sigma_ypE2HE3int_G2_SOS_RasD, 2);
            break;
        case 506:
            dJydy[0] = (-1.0*mypHE3Ev3int_G2_SOS_RasD + 1.0*ypHE3Ev3int_G2_SOS_RasD)/std::pow(sigma_ypHE3Ev3int_G2_SOS_RasD, 2);
            break;
        case 507:
            dJydy[0] = (-1.0*mypE1HE3int_G2_SOS_RasD + 1.0*ypE1HE3int_G2_SOS_RasD)/std::pow(sigma_ypE1HE3int_G2_SOS_RasD, 2);
            break;
        case 508:
            dJydy[0] = (-1.0*mypHE3E4int_G2_SOS_RasD + 1.0*ypHE3E4int_G2_SOS_RasD)/std::pow(sigma_ypHE3E4int_G2_SOS_RasD, 2);
            break;
        case 509:
            dJydy[0] = (-1.0*mypHE3HE4int_G2_SOS_RasD + 1.0*ypHE3HE4int_G2_SOS_RasD)/std::pow(sigma_ypHE3HE4int_G2_SOS_RasD, 2);
            break;
        case 510:
            dJydy[0] = (-1.0*mypE2HE4int_G2_SOS_RasD + 1.0*ypE2HE4int_G2_SOS_RasD)/std::pow(sigma_ypE2HE4int_G2_SOS_RasD, 2);
            break;
        case 511:
            dJydy[0] = (-1.0*mypHE4Ev3int_G2_SOS_RasD + 1.0*ypHE4Ev3int_G2_SOS_RasD)/std::pow(sigma_ypHE4Ev3int_G2_SOS_RasD, 2);
            break;
        case 512:
            dJydy[0] = (-1.0*mypE1HE4int_G2_SOS_RasD + 1.0*ypE1HE4int_G2_SOS_RasD)/std::pow(sigma_ypE1HE4int_G2_SOS_RasD, 2);
            break;
        case 513:
            dJydy[0] = (-1.0*mypE3HE4int_G2_SOS_RasD + 1.0*ypE3HE4int_G2_SOS_RasD)/std::pow(sigma_ypE3HE4int_G2_SOS_RasD, 2);
            break;
        case 514:
            dJydy[0] = (-1.0*mypHE4E4int_G2_SOS_RasD + 1.0*ypHE4E4int_G2_SOS_RasD)/std::pow(sigma_ypHE4E4int_G2_SOS_RasD, 2);
            break;
        case 515:
            dJydy[0] = (-1.0*mypHE4HE4int_G2_SOS_RasD + 1.0*ypHE4HE4int_G2_SOS_RasD)/std::pow(sigma_ypHE4HE4int_G2_SOS_RasD, 2);
            break;
        case 516:
            dJydy[0] = (-1.0*mypHGF_Met_Metint_G2_SOS_RasD + 1.0*ypHGF_Met_Metint_G2_SOS_RasD)/std::pow(sigma_ypHGF_Met_Metint_G2_SOS_RasD, 2);
            break;
        case 517:
            dJydy[0] = (-1.0*mypHGF_Met_HGF_Metint_G2_SOS_RasD + 1.0*ypHGF_Met_HGF_Metint_G2_SOS_RasD)/std::pow(sigma_ypHGF_Met_HGF_Metint_G2_SOS_RasD, 2);
            break;
        case 518:
            dJydy[0] = (-1.0*mypPPrPPrint_G2_SOS_RasD + 1.0*ypPPrPPrint_G2_SOS_RasD)/std::pow(sigma_ypPPrPPrint_G2_SOS_RasD, 2);
            break;
        case 519:
            dJydy[0] = (-1.0*mypPPrPrint_G2_SOS_RasD + 1.0*ypPPrPrint_G2_SOS_RasD)/std::pow(sigma_ypPPrPrint_G2_SOS_RasD, 2);
            break;
        case 520:
            dJydy[0] = (-1.0*mypFFrFFrint_G2_SOS_RasD + 1.0*ypFFrFFrint_G2_SOS_RasD)/std::pow(sigma_ypFFrFFrint_G2_SOS_RasD, 2);
            break;
        case 521:
            dJydy[0] = (-1.0*mypFFrFrint_G2_SOS_RasD + 1.0*ypFFrFrint_G2_SOS_RasD)/std::pow(sigma_ypFFrFrint_G2_SOS_RasD, 2);
            break;
        case 522:
            dJydy[0] = (-1.0*mypIIrIr_int_IRS_G2_SOS_RasD + 1.0*ypIIrIr_int_IRS_G2_SOS_RasD)/std::pow(sigma_ypIIrIr_int_IRS_G2_SOS_RasD, 2);
            break;
        case 523:
            dJydy[0] = (-1.0*mypINS_Isr_Isr_int_IRS_G2_SOS_RasD + 1.0*ypINS_Isr_Isr_int_IRS_G2_SOS_RasD)/std::pow(sigma_ypINS_Isr_Isr_int_IRS_G2_SOS_RasD, 2);
            break;
        case 524:
            dJydy[0] = (-1.0*mypIIrIrI_int_IRS_G2_SOS_RasD + 1.0*ypIIrIrI_int_IRS_G2_SOS_RasD)/std::pow(sigma_ypIIrIrI_int_IRS_G2_SOS_RasD, 2);
            break;
        case 525:
            dJydy[0] = (-1.0*mypINS_Isr_Isr_INS_int_IRS_G2_SOS_RasD + 1.0*ypINS_Isr_Isr_INS_int_IRS_G2_SOS_RasD)/std::pow(sigma_ypINS_Isr_Isr_INS_int_IRS_G2_SOS_RasD, 2);
            break;
        case 526:
            dJydy[0] = (-1.0*mypEE1E2_G2_SOS_RasD + 1.0*ypEE1E2_G2_SOS_RasD)/std::pow(sigma_ypEE1E2_G2_SOS_RasD, 2);
            break;
        case 527:
            dJydy[0] = (-1.0*mypEE1Ev3_G2_SOS_RasD + 1.0*ypEE1Ev3_G2_SOS_RasD)/std::pow(sigma_ypEE1Ev3_G2_SOS_RasD, 2);
            break;
        case 528:
            dJydy[0] = (-1.0*mypEE1E1_G2_SOS_RasD + 1.0*ypEE1E1_G2_SOS_RasD)/std::pow(sigma_ypEE1E1_G2_SOS_RasD, 2);
            break;
        case 529:
            dJydy[0] = (-1.0*mypEE1EE1_G2_SOS_RasD + 1.0*ypEE1EE1_G2_SOS_RasD)/std::pow(sigma_ypEE1EE1_G2_SOS_RasD, 2);
            break;
        case 530:
            dJydy[0] = (-1.0*mypEE1E3_G2_SOS_RasD + 1.0*ypEE1E3_G2_SOS_RasD)/std::pow(sigma_ypEE1E3_G2_SOS_RasD, 2);
            break;
        case 531:
            dJydy[0] = (-1.0*mypEE1HE3_G2_SOS_RasD + 1.0*ypEE1HE3_G2_SOS_RasD)/std::pow(sigma_ypEE1HE3_G2_SOS_RasD, 2);
            break;
        case 532:
            dJydy[0] = (-1.0*mypEE1E4_G2_SOS_RasD + 1.0*ypEE1E4_G2_SOS_RasD)/std::pow(sigma_ypEE1E4_G2_SOS_RasD, 2);
            break;
        case 533:
            dJydy[0] = (-1.0*mypEE1HE4_G2_SOS_RasD + 1.0*ypEE1HE4_G2_SOS_RasD)/std::pow(sigma_ypEE1HE4_G2_SOS_RasD, 2);
            break;
        case 534:
            dJydy[0] = (-1.0*mypE2HE3_G2_SOS_RasD + 1.0*ypE2HE3_G2_SOS_RasD)/std::pow(sigma_ypE2HE3_G2_SOS_RasD, 2);
            break;
        case 535:
            dJydy[0] = (-1.0*mypHE3Ev3_G2_SOS_RasD + 1.0*ypHE3Ev3_G2_SOS_RasD)/std::pow(sigma_ypHE3Ev3_G2_SOS_RasD, 2);
            break;
        case 536:
            dJydy[0] = (-1.0*mypE1HE3_G2_SOS_RasD + 1.0*ypE1HE3_G2_SOS_RasD)/std::pow(sigma_ypE1HE3_G2_SOS_RasD, 2);
            break;
        case 537:
            dJydy[0] = (-1.0*mypHE3E4_G2_SOS_RasD + 1.0*ypHE3E4_G2_SOS_RasD)/std::pow(sigma_ypHE3E4_G2_SOS_RasD, 2);
            break;
        case 538:
            dJydy[0] = (-1.0*mypHE3HE4_G2_SOS_RasD + 1.0*ypHE3HE4_G2_SOS_RasD)/std::pow(sigma_ypHE3HE4_G2_SOS_RasD, 2);
            break;
        case 539:
            dJydy[0] = (-1.0*mypE2HE4_G2_SOS_RasD + 1.0*ypE2HE4_G2_SOS_RasD)/std::pow(sigma_ypE2HE4_G2_SOS_RasD, 2);
            break;
        case 540:
            dJydy[0] = (-1.0*mypHE4Ev3_G2_SOS_RasD + 1.0*ypHE4Ev3_G2_SOS_RasD)/std::pow(sigma_ypHE4Ev3_G2_SOS_RasD, 2);
            break;
        case 541:
            dJydy[0] = (-1.0*mypE1HE4_G2_SOS_RasD + 1.0*ypE1HE4_G2_SOS_RasD)/std::pow(sigma_ypE1HE4_G2_SOS_RasD, 2);
            break;
        case 542:
            dJydy[0] = (-1.0*mypE3HE4_G2_SOS_RasD + 1.0*ypE3HE4_G2_SOS_RasD)/std::pow(sigma_ypE3HE4_G2_SOS_RasD, 2);
            break;
        case 543:
            dJydy[0] = (-1.0*mypHE4E4_G2_SOS_RasD + 1.0*ypHE4E4_G2_SOS_RasD)/std::pow(sigma_ypHE4E4_G2_SOS_RasD, 2);
            break;
        case 544:
            dJydy[0] = (-1.0*mypHE4HE4_G2_SOS_RasD + 1.0*ypHE4HE4_G2_SOS_RasD)/std::pow(sigma_ypHE4HE4_G2_SOS_RasD, 2);
            break;
        case 545:
            dJydy[0] = (-1.0*mypHGF_Met_Met_G2_SOS_RasD + 1.0*ypHGF_Met_Met_G2_SOS_RasD)/std::pow(sigma_ypHGF_Met_Met_G2_SOS_RasD, 2);
            break;
        case 546:
            dJydy[0] = (-1.0*mypHGF_Met_HGF_Met_G2_SOS_RasD + 1.0*ypHGF_Met_HGF_Met_G2_SOS_RasD)/std::pow(sigma_ypHGF_Met_HGF_Met_G2_SOS_RasD, 2);
            break;
        case 547:
            dJydy[0] = (-1.0*mypPPrPPr_G2_SOS_RasD + 1.0*ypPPrPPr_G2_SOS_RasD)/std::pow(sigma_ypPPrPPr_G2_SOS_RasD, 2);
            break;
        case 548:
            dJydy[0] = (-1.0*mypPPrPr_G2_SOS_RasD + 1.0*ypPPrPr_G2_SOS_RasD)/std::pow(sigma_ypPPrPr_G2_SOS_RasD, 2);
            break;
        case 549:
            dJydy[0] = (-1.0*mypFFrFFr_G2_SOS_RasD + 1.0*ypFFrFFr_G2_SOS_RasD)/std::pow(sigma_ypFFrFFr_G2_SOS_RasD, 2);
            break;
        case 550:
            dJydy[0] = (-1.0*mypFFrFr_G2_SOS_RasD + 1.0*ypFFrFr_G2_SOS_RasD)/std::pow(sigma_ypFFrFr_G2_SOS_RasD, 2);
            break;
        case 551:
            dJydy[0] = (-1.0*mypIIrIr_IRS_G2_SOS_RasD + 1.0*ypIIrIr_IRS_G2_SOS_RasD)/std::pow(sigma_ypIIrIr_IRS_G2_SOS_RasD, 2);
            break;
        case 552:
            dJydy[0] = (-1.0*mypINS_Isr_Isr_IRS_G2_SOS_RasD + 1.0*ypINS_Isr_Isr_IRS_G2_SOS_RasD)/std::pow(sigma_ypINS_Isr_Isr_IRS_G2_SOS_RasD, 2);
            break;
        case 553:
            dJydy[0] = (-1.0*mypIIrIrI_IRS_G2_SOS_RasD + 1.0*ypIIrIrI_IRS_G2_SOS_RasD)/std::pow(sigma_ypIIrIrI_IRS_G2_SOS_RasD, 2);
            break;
        case 554:
            dJydy[0] = (-1.0*mypINS_Isr_Isr_INS_IRS_G2_SOS_RasD + 1.0*ypINS_Isr_Isr_INS_IRS_G2_SOS_RasD)/std::pow(sigma_ypINS_Isr_Isr_INS_IRS_G2_SOS_RasD, 2);
            break;
        case 555:
            dJydy[0] = (-1.0*mypEE1E2_PLCg_PIP2 + 1.0*ypEE1E2_PLCg_PIP2)/std::pow(sigma_ypEE1E2_PLCg_PIP2, 2);
            break;
        case 556:
            dJydy[0] = (-1.0*mypEE1Ev3_PLCg_PIP2 + 1.0*ypEE1Ev3_PLCg_PIP2)/std::pow(sigma_ypEE1Ev3_PLCg_PIP2, 2);
            break;
        case 557:
            dJydy[0] = (-1.0*mypEE1E1_PLCg_PIP2 + 1.0*ypEE1E1_PLCg_PIP2)/std::pow(sigma_ypEE1E1_PLCg_PIP2, 2);
            break;
        case 558:
            dJydy[0] = (-1.0*mypEE1EE1_PLCg_PIP2 + 1.0*ypEE1EE1_PLCg_PIP2)/std::pow(sigma_ypEE1EE1_PLCg_PIP2, 2);
            break;
        case 559:
            dJydy[0] = (-1.0*mypEE1E3_PLCg_PIP2 + 1.0*ypEE1E3_PLCg_PIP2)/std::pow(sigma_ypEE1E3_PLCg_PIP2, 2);
            break;
        case 560:
            dJydy[0] = (-1.0*mypEE1HE3_PLCg_PIP2 + 1.0*ypEE1HE3_PLCg_PIP2)/std::pow(sigma_ypEE1HE3_PLCg_PIP2, 2);
            break;
        case 561:
            dJydy[0] = (-1.0*mypEE1E4_PLCg_PIP2 + 1.0*ypEE1E4_PLCg_PIP2)/std::pow(sigma_ypEE1E4_PLCg_PIP2, 2);
            break;
        case 562:
            dJydy[0] = (-1.0*mypEE1HE4_PLCg_PIP2 + 1.0*ypEE1HE4_PLCg_PIP2)/std::pow(sigma_ypEE1HE4_PLCg_PIP2, 2);
            break;
        case 563:
            dJydy[0] = (-1.0*mypE2HE3_PLCg_PIP2 + 1.0*ypE2HE3_PLCg_PIP2)/std::pow(sigma_ypE2HE3_PLCg_PIP2, 2);
            break;
        case 564:
            dJydy[0] = (-1.0*mypHE3Ev3_PLCg_PIP2 + 1.0*ypHE3Ev3_PLCg_PIP2)/std::pow(sigma_ypHE3Ev3_PLCg_PIP2, 2);
            break;
        case 565:
            dJydy[0] = (-1.0*mypE1HE3_PLCg_PIP2 + 1.0*ypE1HE3_PLCg_PIP2)/std::pow(sigma_ypE1HE3_PLCg_PIP2, 2);
            break;
        case 566:
            dJydy[0] = (-1.0*mypHE3E4_PLCg_PIP2 + 1.0*ypHE3E4_PLCg_PIP2)/std::pow(sigma_ypHE3E4_PLCg_PIP2, 2);
            break;
        case 567:
            dJydy[0] = (-1.0*mypHE3HE4_PLCg_PIP2 + 1.0*ypHE3HE4_PLCg_PIP2)/std::pow(sigma_ypHE3HE4_PLCg_PIP2, 2);
            break;
        case 568:
            dJydy[0] = (-1.0*mypE2HE4_PLCg_PIP2 + 1.0*ypE2HE4_PLCg_PIP2)/std::pow(sigma_ypE2HE4_PLCg_PIP2, 2);
            break;
        case 569:
            dJydy[0] = (-1.0*mypHE4Ev3_PLCg_PIP2 + 1.0*ypHE4Ev3_PLCg_PIP2)/std::pow(sigma_ypHE4Ev3_PLCg_PIP2, 2);
            break;
        case 570:
            dJydy[0] = (-1.0*mypE1HE4_PLCg_PIP2 + 1.0*ypE1HE4_PLCg_PIP2)/std::pow(sigma_ypE1HE4_PLCg_PIP2, 2);
            break;
        case 571:
            dJydy[0] = (-1.0*mypE3HE4_PLCg_PIP2 + 1.0*ypE3HE4_PLCg_PIP2)/std::pow(sigma_ypE3HE4_PLCg_PIP2, 2);
            break;
        case 572:
            dJydy[0] = (-1.0*mypHE4E4_PLCg_PIP2 + 1.0*ypHE4E4_PLCg_PIP2)/std::pow(sigma_ypHE4E4_PLCg_PIP2, 2);
            break;
        case 573:
            dJydy[0] = (-1.0*mypHE4HE4_PLCg_PIP2 + 1.0*ypHE4HE4_PLCg_PIP2)/std::pow(sigma_ypHE4HE4_PLCg_PIP2, 2);
            break;
        case 574:
            dJydy[0] = (-1.0*mypHGF_Met_Met_PLCg_PIP2 + 1.0*ypHGF_Met_Met_PLCg_PIP2)/std::pow(sigma_ypHGF_Met_Met_PLCg_PIP2, 2);
            break;
        case 575:
            dJydy[0] = (-1.0*mypHGF_Met_HGF_Met_PLCg_PIP2 + 1.0*ypHGF_Met_HGF_Met_PLCg_PIP2)/std::pow(sigma_ypHGF_Met_HGF_Met_PLCg_PIP2, 2);
            break;
        case 576:
            dJydy[0] = (-1.0*mypPPrPPr_PLCg_PIP2 + 1.0*ypPPrPPr_PLCg_PIP2)/std::pow(sigma_ypPPrPPr_PLCg_PIP2, 2);
            break;
        case 577:
            dJydy[0] = (-1.0*mypPPrPr_PLCg_PIP2 + 1.0*ypPPrPr_PLCg_PIP2)/std::pow(sigma_ypPPrPr_PLCg_PIP2, 2);
            break;
        case 578:
            dJydy[0] = (-1.0*mypFFrFFr_PLCg_PIP2 + 1.0*ypFFrFFr_PLCg_PIP2)/std::pow(sigma_ypFFrFFr_PLCg_PIP2, 2);
            break;
        case 579:
            dJydy[0] = (-1.0*mypFFrFr_PLCg_PIP2 + 1.0*ypFFrFr_PLCg_PIP2)/std::pow(sigma_ypFFrFr_PLCg_PIP2, 2);
            break;
        case 580:
            dJydy[0] = (-1.0*mypIIrIr_IRS_PLCg_PIP2 + 1.0*ypIIrIr_IRS_PLCg_PIP2)/std::pow(sigma_ypIIrIr_IRS_PLCg_PIP2, 2);
            break;
        case 581:
            dJydy[0] = (-1.0*mypINS_Isr_Isr_IRS_PLCg_PIP2 + 1.0*ypINS_Isr_Isr_IRS_PLCg_PIP2)/std::pow(sigma_ypINS_Isr_Isr_IRS_PLCg_PIP2, 2);
            break;
        case 582:
            dJydy[0] = (-1.0*mypIIrIrI_IRS_PLCg_PIP2 + 1.0*ypIIrIrI_IRS_PLCg_PIP2)/std::pow(sigma_ypIIrIrI_IRS_PLCg_PIP2, 2);
            break;
        case 583:
            dJydy[0] = (-1.0*mypINS_Isr_Isr_INS_IRS_PLCg_PIP2 + 1.0*ypINS_Isr_Isr_INS_IRS_PLCg_PIP2)/std::pow(sigma_ypINS_Isr_Isr_INS_IRS_PLCg_PIP2, 2);
            break;
        case 584:
            dJydy[0] = (-1.0*mypEE1E2_PI3K1_PIP2 + 1.0*ypEE1E2_PI3K1_PIP2)/std::pow(sigma_ypEE1E2_PI3K1_PIP2, 2);
            break;
        case 585:
            dJydy[0] = (-1.0*mypEE1Ev3_PI3K1_PIP2 + 1.0*ypEE1Ev3_PI3K1_PIP2)/std::pow(sigma_ypEE1Ev3_PI3K1_PIP2, 2);
            break;
        case 586:
            dJydy[0] = (-1.0*mypEE1E1_PI3K1_PIP2 + 1.0*ypEE1E1_PI3K1_PIP2)/std::pow(sigma_ypEE1E1_PI3K1_PIP2, 2);
            break;
        case 587:
            dJydy[0] = (-1.0*mypEE1EE1_PI3K1_PIP2 + 1.0*ypEE1EE1_PI3K1_PIP2)/std::pow(sigma_ypEE1EE1_PI3K1_PIP2, 2);
            break;
        case 588:
            dJydy[0] = (-1.0*mypEE1E3_PI3K1_PIP2 + 1.0*ypEE1E3_PI3K1_PIP2)/std::pow(sigma_ypEE1E3_PI3K1_PIP2, 2);
            break;
        case 589:
            dJydy[0] = (-1.0*mypEE1HE3_PI3K1_PIP2 + 1.0*ypEE1HE3_PI3K1_PIP2)/std::pow(sigma_ypEE1HE3_PI3K1_PIP2, 2);
            break;
        case 590:
            dJydy[0] = (-1.0*mypEE1E4_PI3K1_PIP2 + 1.0*ypEE1E4_PI3K1_PIP2)/std::pow(sigma_ypEE1E4_PI3K1_PIP2, 2);
            break;
        case 591:
            dJydy[0] = (-1.0*mypEE1HE4_PI3K1_PIP2 + 1.0*ypEE1HE4_PI3K1_PIP2)/std::pow(sigma_ypEE1HE4_PI3K1_PIP2, 2);
            break;
        case 592:
            dJydy[0] = (-1.0*mypE2HE3_PI3K1_PIP2 + 1.0*ypE2HE3_PI3K1_PIP2)/std::pow(sigma_ypE2HE3_PI3K1_PIP2, 2);
            break;
        case 593:
            dJydy[0] = (-1.0*mypHE3Ev3_PI3K1_PIP2 + 1.0*ypHE3Ev3_PI3K1_PIP2)/std::pow(sigma_ypHE3Ev3_PI3K1_PIP2, 2);
            break;
        case 594:
            dJydy[0] = (-1.0*mypE1HE3_PI3K1_PIP2 + 1.0*ypE1HE3_PI3K1_PIP2)/std::pow(sigma_ypE1HE3_PI3K1_PIP2, 2);
            break;
        case 595:
            dJydy[0] = (-1.0*mypHE3E4_PI3K1_PIP2 + 1.0*ypHE3E4_PI3K1_PIP2)/std::pow(sigma_ypHE3E4_PI3K1_PIP2, 2);
            break;
        case 596:
            dJydy[0] = (-1.0*mypHE3HE4_PI3K1_PIP2 + 1.0*ypHE3HE4_PI3K1_PIP2)/std::pow(sigma_ypHE3HE4_PI3K1_PIP2, 2);
            break;
        case 597:
            dJydy[0] = (-1.0*mypE2HE4_PI3K1_PIP2 + 1.0*ypE2HE4_PI3K1_PIP2)/std::pow(sigma_ypE2HE4_PI3K1_PIP2, 2);
            break;
        case 598:
            dJydy[0] = (-1.0*mypHE4Ev3_PI3K1_PIP2 + 1.0*ypHE4Ev3_PI3K1_PIP2)/std::pow(sigma_ypHE4Ev3_PI3K1_PIP2, 2);
            break;
        case 599:
            dJydy[0] = (-1.0*mypE1HE4_PI3K1_PIP2 + 1.0*ypE1HE4_PI3K1_PIP2)/std::pow(sigma_ypE1HE4_PI3K1_PIP2, 2);
            break;
        case 600:
            dJydy[0] = (-1.0*mypE3HE4_PI3K1_PIP2 + 1.0*ypE3HE4_PI3K1_PIP2)/std::pow(sigma_ypE3HE4_PI3K1_PIP2, 2);
            break;
        case 601:
            dJydy[0] = (-1.0*mypHE4E4_PI3K1_PIP2 + 1.0*ypHE4E4_PI3K1_PIP2)/std::pow(sigma_ypHE4E4_PI3K1_PIP2, 2);
            break;
        case 602:
            dJydy[0] = (-1.0*mypHE4HE4_PI3K1_PIP2 + 1.0*ypHE4HE4_PI3K1_PIP2)/std::pow(sigma_ypHE4HE4_PI3K1_PIP2, 2);
            break;
        case 603:
            dJydy[0] = (-1.0*mypHGF_Met_Met_PI3K1_PIP2 + 1.0*ypHGF_Met_Met_PI3K1_PIP2)/std::pow(sigma_ypHGF_Met_Met_PI3K1_PIP2, 2);
            break;
        case 604:
            dJydy[0] = (-1.0*mypHGF_Met_HGF_Met_PI3K1_PIP2 + 1.0*ypHGF_Met_HGF_Met_PI3K1_PIP2)/std::pow(sigma_ypHGF_Met_HGF_Met_PI3K1_PIP2, 2);
            break;
        case 605:
            dJydy[0] = (-1.0*mypPPrPPr_PI3K1_PIP2 + 1.0*ypPPrPPr_PI3K1_PIP2)/std::pow(sigma_ypPPrPPr_PI3K1_PIP2, 2);
            break;
        case 606:
            dJydy[0] = (-1.0*mypPPrPr_PI3K1_PIP2 + 1.0*ypPPrPr_PI3K1_PIP2)/std::pow(sigma_ypPPrPr_PI3K1_PIP2, 2);
            break;
        case 607:
            dJydy[0] = (-1.0*mypFFrFFr_PI3K1_PIP2 + 1.0*ypFFrFFr_PI3K1_PIP2)/std::pow(sigma_ypFFrFFr_PI3K1_PIP2, 2);
            break;
        case 608:
            dJydy[0] = (-1.0*mypFFrFr_PI3K1_PIP2 + 1.0*ypFFrFr_PI3K1_PIP2)/std::pow(sigma_ypFFrFr_PI3K1_PIP2, 2);
            break;
        case 609:
            dJydy[0] = (-1.0*mypIIrIr_IRS_PI3K1_PIP2 + 1.0*ypIIrIr_IRS_PI3K1_PIP2)/std::pow(sigma_ypIIrIr_IRS_PI3K1_PIP2, 2);
            break;
        case 610:
            dJydy[0] = (-1.0*mypINS_Isr_Isr_IRS_PI3K1_PIP2 + 1.0*ypINS_Isr_Isr_IRS_PI3K1_PIP2)/std::pow(sigma_ypINS_Isr_Isr_IRS_PI3K1_PIP2, 2);
            break;
        case 611:
            dJydy[0] = (-1.0*mypIIrIrI_IRS_PI3K1_PIP2 + 1.0*ypIIrIrI_IRS_PI3K1_PIP2)/std::pow(sigma_ypIIrIrI_IRS_PI3K1_PIP2, 2);
            break;
        case 612:
            dJydy[0] = (-1.0*mypINS_Isr_Isr_INS_IRS_PI3K1_PIP2 + 1.0*ypINS_Isr_Isr_INS_IRS_PI3K1_PIP2)/std::pow(sigma_ypINS_Isr_Isr_INS_IRS_PI3K1_PIP2, 2);
            break;
        case 613:
            dJydy[0] = (-1.0*mypEE1E2_PI3K2_PIP + 1.0*ypEE1E2_PI3K2_PIP)/std::pow(sigma_ypEE1E2_PI3K2_PIP, 2);
            break;
        case 614:
            dJydy[0] = (-1.0*mypEE1Ev3_PI3K2_PIP + 1.0*ypEE1Ev3_PI3K2_PIP)/std::pow(sigma_ypEE1Ev3_PI3K2_PIP, 2);
            break;
        case 615:
            dJydy[0] = (-1.0*mypEE1E1_PI3K2_PIP + 1.0*ypEE1E1_PI3K2_PIP)/std::pow(sigma_ypEE1E1_PI3K2_PIP, 2);
            break;
        case 616:
            dJydy[0] = (-1.0*mypEE1EE1_PI3K2_PIP + 1.0*ypEE1EE1_PI3K2_PIP)/std::pow(sigma_ypEE1EE1_PI3K2_PIP, 2);
            break;
        case 617:
            dJydy[0] = (-1.0*mypEE1E3_PI3K2_PIP + 1.0*ypEE1E3_PI3K2_PIP)/std::pow(sigma_ypEE1E3_PI3K2_PIP, 2);
            break;
        case 618:
            dJydy[0] = (-1.0*mypEE1HE3_PI3K2_PIP + 1.0*ypEE1HE3_PI3K2_PIP)/std::pow(sigma_ypEE1HE3_PI3K2_PIP, 2);
            break;
        case 619:
            dJydy[0] = (-1.0*mypEE1E4_PI3K2_PIP + 1.0*ypEE1E4_PI3K2_PIP)/std::pow(sigma_ypEE1E4_PI3K2_PIP, 2);
            break;
        case 620:
            dJydy[0] = (-1.0*mypEE1HE4_PI3K2_PIP + 1.0*ypEE1HE4_PI3K2_PIP)/std::pow(sigma_ypEE1HE4_PI3K2_PIP, 2);
            break;
        case 621:
            dJydy[0] = (-1.0*mypE2HE3_PI3K2_PIP + 1.0*ypE2HE3_PI3K2_PIP)/std::pow(sigma_ypE2HE3_PI3K2_PIP, 2);
            break;
        case 622:
            dJydy[0] = (-1.0*mypHE3Ev3_PI3K2_PIP + 1.0*ypHE3Ev3_PI3K2_PIP)/std::pow(sigma_ypHE3Ev3_PI3K2_PIP, 2);
            break;
        case 623:
            dJydy[0] = (-1.0*mypE1HE3_PI3K2_PIP + 1.0*ypE1HE3_PI3K2_PIP)/std::pow(sigma_ypE1HE3_PI3K2_PIP, 2);
            break;
        case 624:
            dJydy[0] = (-1.0*mypHE3E4_PI3K2_PIP + 1.0*ypHE3E4_PI3K2_PIP)/std::pow(sigma_ypHE3E4_PI3K2_PIP, 2);
            break;
        case 625:
            dJydy[0] = (-1.0*mypHE3HE4_PI3K2_PIP + 1.0*ypHE3HE4_PI3K2_PIP)/std::pow(sigma_ypHE3HE4_PI3K2_PIP, 2);
            break;
        case 626:
            dJydy[0] = (-1.0*mypE2HE4_PI3K2_PIP + 1.0*ypE2HE4_PI3K2_PIP)/std::pow(sigma_ypE2HE4_PI3K2_PIP, 2);
            break;
        case 627:
            dJydy[0] = (-1.0*mypHE4Ev3_PI3K2_PIP + 1.0*ypHE4Ev3_PI3K2_PIP)/std::pow(sigma_ypHE4Ev3_PI3K2_PIP, 2);
            break;
        case 628:
            dJydy[0] = (-1.0*mypE1HE4_PI3K2_PIP + 1.0*ypE1HE4_PI3K2_PIP)/std::pow(sigma_ypE1HE4_PI3K2_PIP, 2);
            break;
        case 629:
            dJydy[0] = (-1.0*mypE3HE4_PI3K2_PIP + 1.0*ypE3HE4_PI3K2_PIP)/std::pow(sigma_ypE3HE4_PI3K2_PIP, 2);
            break;
        case 630:
            dJydy[0] = (-1.0*mypHE4E4_PI3K2_PIP + 1.0*ypHE4E4_PI3K2_PIP)/std::pow(sigma_ypHE4E4_PI3K2_PIP, 2);
            break;
        case 631:
            dJydy[0] = (-1.0*mypHE4HE4_PI3K2_PIP + 1.0*ypHE4HE4_PI3K2_PIP)/std::pow(sigma_ypHE4HE4_PI3K2_PIP, 2);
            break;
        case 632:
            dJydy[0] = (-1.0*mypHGF_Met_Met_PI3K2_PIP + 1.0*ypHGF_Met_Met_PI3K2_PIP)/std::pow(sigma_ypHGF_Met_Met_PI3K2_PIP, 2);
            break;
        case 633:
            dJydy[0] = (-1.0*mypHGF_Met_HGF_Met_PI3K2_PIP + 1.0*ypHGF_Met_HGF_Met_PI3K2_PIP)/std::pow(sigma_ypHGF_Met_HGF_Met_PI3K2_PIP, 2);
            break;
        case 634:
            dJydy[0] = (-1.0*mypPPrPPr_PI3K2_PIP + 1.0*ypPPrPPr_PI3K2_PIP)/std::pow(sigma_ypPPrPPr_PI3K2_PIP, 2);
            break;
        case 635:
            dJydy[0] = (-1.0*mypPPrPr_PI3K2_PIP + 1.0*ypPPrPr_PI3K2_PIP)/std::pow(sigma_ypPPrPr_PI3K2_PIP, 2);
            break;
        case 636:
            dJydy[0] = (-1.0*mypFFrFFr_PI3K2_PIP + 1.0*ypFFrFFr_PI3K2_PIP)/std::pow(sigma_ypFFrFFr_PI3K2_PIP, 2);
            break;
        case 637:
            dJydy[0] = (-1.0*mypFFrFr_PI3K2_PIP + 1.0*ypFFrFr_PI3K2_PIP)/std::pow(sigma_ypFFrFr_PI3K2_PIP, 2);
            break;
        case 638:
            dJydy[0] = (-1.0*mypIIrIr_IRS_PI3K2_PIP + 1.0*ypIIrIr_IRS_PI3K2_PIP)/std::pow(sigma_ypIIrIr_IRS_PI3K2_PIP, 2);
            break;
        case 639:
            dJydy[0] = (-1.0*mypINS_Isr_Isr_IRS_PI3K2_PIP + 1.0*ypINS_Isr_Isr_IRS_PI3K2_PIP)/std::pow(sigma_ypINS_Isr_Isr_IRS_PI3K2_PIP, 2);
            break;
        case 640:
            dJydy[0] = (-1.0*mypIIrIrI_IRS_PI3K2_PIP + 1.0*ypIIrIrI_IRS_PI3K2_PIP)/std::pow(sigma_ypIIrIrI_IRS_PI3K2_PIP, 2);
            break;
        case 641:
            dJydy[0] = (-1.0*mypINS_Isr_Isr_INS_IRS_PI3K2_PIP + 1.0*ypINS_Isr_Isr_INS_IRS_PI3K2_PIP)/std::pow(sigma_ypINS_Isr_Isr_INS_IRS_PI3K2_PIP, 2);
            break;
        case 642:
            dJydy[0] = (-1.0*myIRS + 1.0*yIRS)/std::pow(sigma_yIRS, 2);
            break;
        case 643:
            dJydy[0] = (-1.0*mySp + 1.0*ySp)/std::pow(sigma_ySp, 2);
            break;
        case 644:
            dJydy[0] = (-1.0*myCbl + 1.0*yCbl)/std::pow(sigma_yCbl, 2);
            break;
        case 645:
            dJydy[0] = (-1.0*myG2 + 1.0*yG2)/std::pow(sigma_yG2, 2);
            break;
        case 646:
            dJydy[0] = (-1.0*myG2_SOS + 1.0*yG2_SOS)/std::pow(sigma_yG2_SOS, 2);
            break;
        case 647:
            dJydy[0] = (-1.0*myG2_pSOS + 1.0*yG2_pSOS)/std::pow(sigma_yG2_pSOS, 2);
            break;
        case 648:
            dJydy[0] = (-1.0*myPLCg + 1.0*yPLCg)/std::pow(sigma_yPLCg, 2);
            break;
        case 649:
            dJydy[0] = (-1.0*myPI3KC1 + 1.0*yPI3KC1)/std::pow(sigma_yPI3KC1, 2);
            break;
        case 650:
            dJydy[0] = (-1.0*myPI3KR1 + 1.0*yPI3KR1)/std::pow(sigma_yPI3KR1, 2);
            break;
        case 651:
            dJydy[0] = (-1.0*myPI3K1 + 1.0*yPI3K1)/std::pow(sigma_yPI3K1, 2);
            break;
        case 652:
            dJydy[0] = (-1.0*mypPI3K1 + 1.0*ypPI3K1)/std::pow(sigma_ypPI3K1, 2);
            break;
        case 653:
            dJydy[0] = (-1.0*myPI3K2 + 1.0*yPI3K2)/std::pow(sigma_yPI3K2, 2);
            break;
        case 654:
            dJydy[0] = (-1.0*mymTORC1 + 1.0*ymTORC1)/std::pow(sigma_ymTORC1, 2);
            break;
        case 655:
            dJydy[0] = (-1.0*mymTORC1active + 1.0*ymTORC1active)/std::pow(sigma_ymTORC1active, 2);
            break;
        case 656:
            dJydy[0] = (-1.0*myPIP + 1.0*yPIP)/std::pow(sigma_yPIP, 2);
            break;
        case 657:
            dJydy[0] = (-1.0*myPI3P + 1.0*yPI3P)/std::pow(sigma_yPI3P, 2);
            break;
        case 658:
            dJydy[0] = (-1.0*myDAG + 1.0*yDAG)/std::pow(sigma_yDAG, 2);
            break;
        case 659:
            dJydy[0] = (-1.0*myGRP + 1.0*yGRP)/std::pow(sigma_yGRP, 2);
            break;
        case 660:
            dJydy[0] = (-1.0*myDAG_GRP + 1.0*yDAG_GRP)/std::pow(sigma_yDAG_GRP, 2);
            break;
        case 661:
            dJydy[0] = (-1.0*myRasT + 1.0*yRasT)/std::pow(sigma_yRasT, 2);
            break;
        case 662:
            dJydy[0] = (-1.0*myRasD + 1.0*yRasD)/std::pow(sigma_yRasD, 2);
            break;
        case 663:
            dJydy[0] = (-1.0*myNF1 + 1.0*yNF1)/std::pow(sigma_yNF1, 2);
            break;
        case 664:
            dJydy[0] = (-1.0*mypNF1 + 1.0*ypNF1)/std::pow(sigma_ypNF1, 2);
            break;
        case 665:
            dJydy[0] = (-1.0*mypCRaf + 1.0*ypCRaf)/std::pow(sigma_ypCRaf, 2);
            break;
        case 666:
            dJydy[0] = (-1.0*myCRaf + 1.0*yCRaf)/std::pow(sigma_yCRaf, 2);
            break;
        case 667:
            dJydy[0] = (-1.0*myRasT_CRaf + 1.0*yRasT_CRaf)/std::pow(sigma_yRasT_CRaf, 2);
            break;
        case 668:
            dJydy[0] = (-1.0*myBRaf + 1.0*yBRaf)/std::pow(sigma_yBRaf, 2);
            break;
        case 669:
            dJydy[0] = (-1.0*myRasT_CRaf_BRaf + 1.0*yRasT_CRaf_BRaf)/std::pow(sigma_yRasT_CRaf_BRaf, 2);
            break;
        case 670:
            dJydy[0] = (-1.0*myMEK + 1.0*yMEK)/std::pow(sigma_yMEK, 2);
            break;
        case 671:
            dJydy[0] = (-1.0*mypMEK + 1.0*ypMEK)/std::pow(sigma_ypMEK, 2);
            break;
        case 672:
            dJydy[0] = (-1.0*myppMEK + 1.0*yppMEK)/std::pow(sigma_yppMEK, 2);
            break;
        case 673:
            dJydy[0] = (-1.0*myMKP3 + 1.0*yMKP3)/std::pow(sigma_yMKP3, 2);
            break;
        case 674:
            dJydy[0] = (-1.0*myERKnuc + 1.0*yERKnuc)/std::pow(sigma_yERKnuc, 2);
            break;
        case 675:
            dJydy[0] = (-1.0*myppERKnuc + 1.0*yppERKnuc)/std::pow(sigma_yppERKnuc, 2);
            break;
        case 676:
            dJydy[0] = (-1.0*myRSK + 1.0*yRSK)/std::pow(sigma_yRSK, 2);
            break;
        case 677:
            dJydy[0] = (-1.0*mypRSK + 1.0*ypRSK)/std::pow(sigma_ypRSK, 2);
            break;
        case 678:
            dJydy[0] = (-1.0*mypRSKnuc + 1.0*ypRSKnuc)/std::pow(sigma_ypRSKnuc, 2);
            break;
        case 679:
            dJydy[0] = (-1.0*myMKP1 + 1.0*yMKP1)/std::pow(sigma_yMKP1, 2);
            break;
        case 680:
            dJydy[0] = (-1.0*mypMKP1 + 1.0*ypMKP1)/std::pow(sigma_ypMKP1, 2);
            break;
        case 681:
            dJydy[0] = (-1.0*mycFos + 1.0*ycFos)/std::pow(sigma_ycFos, 2);
            break;
        case 682:
            dJydy[0] = (-1.0*mypcFos + 1.0*ypcFos)/std::pow(sigma_ypcFos, 2);
            break;
        case 683:
            dJydy[0] = (-1.0*mycJun + 1.0*ycJun)/std::pow(sigma_ycJun, 2);
            break;
        case 684:
            dJydy[0] = (-1.0*mypcFos_cJun + 1.0*ypcFos_cJun)/std::pow(sigma_ypcFos_cJun, 2);
            break;
        case 685:
            dJydy[0] = (-1.0*mycMyc + 1.0*ycMyc)/std::pow(sigma_ycMyc, 2);
            break;
        case 686:
            dJydy[0] = (-1.0*mybCATENINnuc + 1.0*ybCATENINnuc)/std::pow(sigma_ybCATENINnuc, 2);
            break;
        case 687:
            dJydy[0] = (-1.0*mybCATENIN + 1.0*ybCATENIN)/std::pow(sigma_ybCATENIN, 2);
            break;
        case 688:
            dJydy[0] = (-1.0*mypbCATENIN + 1.0*ypbCATENIN)/std::pow(sigma_ypbCATENIN, 2);
            break;
        case 689:
            dJydy[0] = (-1.0*myIP3 + 1.0*yIP3)/std::pow(sigma_yIP3, 2);
            break;
        case 690:
            dJydy[0] = (-1.0*myPIP2 + 1.0*yPIP2)/std::pow(sigma_yPIP2, 2);
            break;
        case 691:
            dJydy[0] = (-1.0*myPIP3 + 1.0*yPIP3)/std::pow(sigma_yPIP3, 2);
            break;
        case 692:
            dJydy[0] = (-1.0*myPTEN + 1.0*yPTEN)/std::pow(sigma_yPTEN, 2);
            break;
        case 693:
            dJydy[0] = (-1.0*myPIP3_AKT + 1.0*yPIP3_AKT)/std::pow(sigma_yPIP3_AKT, 2);
            break;
        case 694:
            dJydy[0] = (-1.0*myAKT + 1.0*yAKT)/std::pow(sigma_yAKT, 2);
            break;
        case 695:
            dJydy[0] = (-1.0*mypAKT + 1.0*ypAKT)/std::pow(sigma_ypAKT, 2);
            break;
        case 696:
            dJydy[0] = (-1.0*myppAKT + 1.0*yppAKT)/std::pow(sigma_yppAKT, 2);
            break;
        case 697:
            dJydy[0] = (-1.0*myPDK1 + 1.0*yPDK1)/std::pow(sigma_yPDK1, 2);
            break;
        case 698:
            dJydy[0] = (-1.0*myPIP3_PDK1 + 1.0*yPIP3_PDK1)/std::pow(sigma_yPIP3_PDK1, 2);
            break;
        case 699:
            dJydy[0] = (-1.0*myPIP3_pAKT + 1.0*yPIP3_pAKT)/std::pow(sigma_yPIP3_pAKT, 2);
            break;
        case 700:
            dJydy[0] = (-1.0*myRictor + 1.0*yRictor)/std::pow(sigma_yRictor, 2);
            break;
        case 701:
            dJydy[0] = (-1.0*mymTOR + 1.0*ymTOR)/std::pow(sigma_ymTOR, 2);
            break;
        case 702:
            dJydy[0] = (-1.0*mymTORC2 + 1.0*ymTORC2)/std::pow(sigma_ymTORC2, 2);
            break;
        case 703:
            dJydy[0] = (-1.0*myPIP3_ppAKT + 1.0*yPIP3_ppAKT)/std::pow(sigma_yPIP3_ppAKT, 2);
            break;
        case 704:
            dJydy[0] = (-1.0*myGSK3b + 1.0*yGSK3b)/std::pow(sigma_yGSK3b, 2);
            break;
        case 705:
            dJydy[0] = (-1.0*mypGSK3b + 1.0*ypGSK3b)/std::pow(sigma_ypGSK3b, 2);
            break;
        case 706:
            dJydy[0] = (-1.0*myTSC1 + 1.0*yTSC1)/std::pow(sigma_yTSC1, 2);
            break;
        case 707:
            dJydy[0] = (-1.0*myTSC2 + 1.0*yTSC2)/std::pow(sigma_yTSC2, 2);
            break;
        case 708:
            dJydy[0] = (-1.0*mypTSC2 + 1.0*ypTSC2)/std::pow(sigma_ypTSC2, 2);
            break;
        case 709:
            dJydy[0] = (-1.0*myTSC + 1.0*yTSC)/std::pow(sigma_yTSC, 2);
            break;
        case 710:
            dJydy[0] = (-1.0*myPKC + 1.0*yPKC)/std::pow(sigma_yPKC, 2);
            break;
        case 711:
            dJydy[0] = (-1.0*myDAG_PKC + 1.0*yDAG_PKC)/std::pow(sigma_yDAG_PKC, 2);
            break;
        case 712:
            dJydy[0] = (-1.0*mypRKIP + 1.0*ypRKIP)/std::pow(sigma_ypRKIP, 2);
            break;
        case 713:
            dJydy[0] = (-1.0*myRKIP + 1.0*yRKIP)/std::pow(sigma_yRKIP, 2);
            break;
        case 714:
            dJydy[0] = (-1.0*myRKIP_CRaf + 1.0*yRKIP_CRaf)/std::pow(sigma_yRKIP_CRaf, 2);
            break;
        case 715:
            dJydy[0] = (-1.0*myERK + 1.0*yERK)/std::pow(sigma_yERK, 2);
            break;
        case 716:
            dJydy[0] = (-1.0*mypERK + 1.0*ypERK)/std::pow(sigma_ypERK, 2);
            break;
        case 717:
            dJydy[0] = (-1.0*myppERK + 1.0*yppERK)/std::pow(sigma_yppERK, 2);
            break;
        case 718:
            dJydy[0] = (-1.0*myFOXO + 1.0*yFOXO)/std::pow(sigma_yFOXO, 2);
            break;
        case 719:
            dJydy[0] = (-1.0*mypFOXO + 1.0*ypFOXO)/std::pow(sigma_ypFOXO, 2);
            break;
        case 720:
            dJydy[0] = (-1.0*myRhebD + 1.0*yRhebD)/std::pow(sigma_yRhebD, 2);
            break;
        case 721:
            dJydy[0] = (-1.0*myRhebT + 1.0*yRhebT)/std::pow(sigma_yRhebT, 2);
            break;
        case 722:
            dJydy[0] = (-1.0*myRaptor + 1.0*yRaptor)/std::pow(sigma_yRaptor, 2);
            break;
        case 723:
            dJydy[0] = (-1.0*myS6K + 1.0*yS6K)/std::pow(sigma_yS6K, 2);
            break;
        case 724:
            dJydy[0] = (-1.0*mypS6K + 1.0*ypS6K)/std::pow(sigma_ypS6K, 2);
            break;
        case 725:
            dJydy[0] = (-1.0*myEIF4EBP1 + 1.0*yEIF4EBP1)/std::pow(sigma_yEIF4EBP1, 2);
            break;
        case 726:
            dJydy[0] = (-1.0*mypEIF4EBP1 + 1.0*ypEIF4EBP1)/std::pow(sigma_ypEIF4EBP1, 2);
            break;
        case 727:
            dJydy[0] = (-1.0*mySOS + 1.0*ySOS)/std::pow(sigma_ySOS, 2);
            break;
        case 728:
            dJydy[0] = (-1.0*myG2_SOS_ppERK + 1.0*yG2_SOS_ppERK)/std::pow(sigma_yG2_SOS_ppERK, 2);
            break;
        case 729:
            dJydy[0] = (-1.0*myCRaf_ppERK + 1.0*yCRaf_ppERK)/std::pow(sigma_yCRaf_ppERK, 2);
            break;
        case 730:
            dJydy[0] = (-1.0*myRasD_DAG_GRP + 1.0*yRasD_DAG_GRP)/std::pow(sigma_yRasD_DAG_GRP, 2);
            break;
        case 731:
            dJydy[0] = (-1.0*myRasT_NF1 + 1.0*yRasT_NF1)/std::pow(sigma_yRasT_NF1, 2);
            break;
        case 732:
            dJydy[0] = (-1.0*myNF1_ppERK + 1.0*yNF1_ppERK)/std::pow(sigma_yNF1_ppERK, 2);
            break;
        case 733:
            dJydy[0] = (-1.0*myMEK_RasT_CRaf_BRaf + 1.0*yMEK_RasT_CRaf_BRaf)/std::pow(sigma_yMEK_RasT_CRaf_BRaf, 2);
            break;
        case 734:
            dJydy[0] = (-1.0*mypMEK_RasT_CRaf_BRaf + 1.0*ypMEK_RasT_CRaf_BRaf)/std::pow(sigma_ypMEK_RasT_CRaf_BRaf, 2);
            break;
        case 735:
            dJydy[0] = (-1.0*myERK_ppMEK + 1.0*yERK_ppMEK)/std::pow(sigma_yERK_ppMEK, 2);
            break;
        case 736:
            dJydy[0] = (-1.0*mypERK_ppMEK + 1.0*ypERK_ppMEK)/std::pow(sigma_ypERK_ppMEK, 2);
            break;
        case 737:
            dJydy[0] = (-1.0*myRSK_ppERK + 1.0*yRSK_ppERK)/std::pow(sigma_yRSK_ppERK, 2);
            break;
        case 738:
            dJydy[0] = (-1.0*mypRSKnuc_MKP1 + 1.0*ypRSKnuc_MKP1)/std::pow(sigma_ypRSKnuc_MKP1, 2);
            break;
        case 739:
            dJydy[0] = (-1.0*myppERKnuc_MKP1 + 1.0*yppERKnuc_MKP1)/std::pow(sigma_yppERKnuc_MKP1, 2);
            break;
        case 740:
            dJydy[0] = (-1.0*mycFos_pRSKnuc + 1.0*ycFos_pRSKnuc)/std::pow(sigma_ycFos_pRSKnuc, 2);
            break;
        case 741:
            dJydy[0] = (-1.0*mycFos_ppERKnuc + 1.0*ycFos_ppERKnuc)/std::pow(sigma_ycFos_ppERKnuc, 2);
            break;
        case 742:
            dJydy[0] = (-1.0*myRKIP_DAG_PKC + 1.0*yRKIP_DAG_PKC)/std::pow(sigma_yRKIP_DAG_PKC, 2);
            break;
        case 743:
            dJydy[0] = (-1.0*myPIP3_PTEN + 1.0*yPIP3_PTEN)/std::pow(sigma_yPIP3_PTEN, 2);
            break;
        case 744:
            dJydy[0] = (-1.0*myPIP3_AKT_PIP3_PDK1 + 1.0*yPIP3_AKT_PIP3_PDK1)/std::pow(sigma_yPIP3_AKT_PIP3_PDK1, 2);
            break;
        case 745:
            dJydy[0] = (-1.0*myPIP3_pAKT_mTORC2 + 1.0*yPIP3_pAKT_mTORC2)/std::pow(sigma_yPIP3_pAKT_mTORC2, 2);
            break;
        case 746:
            dJydy[0] = (-1.0*myGSK3b_ppAKT + 1.0*yGSK3b_ppAKT)/std::pow(sigma_yGSK3b_ppAKT, 2);
            break;
        case 747:
            dJydy[0] = (-1.0*myTSC2_ppAKT + 1.0*yTSC2_ppAKT)/std::pow(sigma_yTSC2_ppAKT, 2);
            break;
        case 748:
            dJydy[0] = (-1.0*myTSC2_ppERK + 1.0*yTSC2_ppERK)/std::pow(sigma_yTSC2_ppERK, 2);
            break;
        case 749:
            dJydy[0] = (-1.0*myRhebT_TSC + 1.0*yRhebT_TSC)/std::pow(sigma_yRhebT_TSC, 2);
            break;
        case 750:
            dJydy[0] = (-1.0*myEIF4EBP1_mTORC1active + 1.0*yEIF4EBP1_mTORC1active)/std::pow(sigma_yEIF4EBP1_mTORC1active, 2);
            break;
        case 751:
            dJydy[0] = (-1.0*myS6K_mTORC1active + 1.0*yS6K_mTORC1active)/std::pow(sigma_yS6K_mTORC1active, 2);
            break;
        case 752:
            dJydy[0] = (-1.0*myFOXO_ppAKT + 1.0*yFOXO_ppAKT)/std::pow(sigma_yFOXO_ppAKT, 2);
            break;
        case 753:
            dJydy[0] = (-1.0*myPI3K1_mTORC1active + 1.0*yPI3K1_mTORC1active)/std::pow(sigma_yPI3K1_mTORC1active, 2);
            break;
        case 754:
            dJydy[0] = (-1.0*mypERK_MKP3 + 1.0*ypERK_MKP3)/std::pow(sigma_ypERK_MKP3, 2);
            break;
        case 755:
            dJydy[0] = (-1.0*myppERK_MKP3 + 1.0*yppERK_MKP3)/std::pow(sigma_yppERK_MKP3, 2);
            break;
        case 756:
            dJydy[0] = (-1.0*myppERKnuc_pMKP1 + 1.0*yppERKnuc_pMKP1)/std::pow(sigma_yppERKnuc_pMKP1, 2);
            break;
        case 757:
            dJydy[0] = (-1.0*myRasT_BRaf + 1.0*yRasT_BRaf)/std::pow(sigma_yRasT_BRaf, 2);
            break;
        case 758:
            dJydy[0] = (-1.0*myRasT_BRaf_BRaf + 1.0*yRasT_BRaf_BRaf)/std::pow(sigma_yRasT_BRaf_BRaf, 2);
            break;
        case 759:
            dJydy[0] = (-1.0*myMEK_RasT_BRaf_BRaf + 1.0*yMEK_RasT_BRaf_BRaf)/std::pow(sigma_yMEK_RasT_BRaf_BRaf, 2);
            break;
        case 760:
            dJydy[0] = (-1.0*mypMEK_RasT_BRaf_BRaf + 1.0*ypMEK_RasT_BRaf_BRaf)/std::pow(sigma_ypMEK_RasT_BRaf_BRaf, 2);
            break;
        case 761:
            dJydy[0] = (-1.0*myEIF4E + 1.0*yEIF4E)/std::pow(sigma_yEIF4E, 2);
            break;
        case 762:
            dJydy[0] = (-1.0*myEIF4EBP1_EIF4E + 1.0*yEIF4EBP1_EIF4E)/std::pow(sigma_yEIF4EBP1_EIF4E, 2);
            break;
        case 763:
            dJydy[0] = (-1.0*myRasT_CRaf_CRaf + 1.0*yRasT_CRaf_CRaf)/std::pow(sigma_yRasT_CRaf_CRaf, 2);
            break;
        case 764:
            dJydy[0] = (-1.0*myMEK_RasT_CRaf_CRaf + 1.0*yMEK_RasT_CRaf_CRaf)/std::pow(sigma_yMEK_RasT_CRaf_CRaf, 2);
            break;
        case 765:
            dJydy[0] = (-1.0*mypMEK_RasT_CRaf_CRaf + 1.0*ypMEK_RasT_CRaf_CRaf)/std::pow(sigma_ypMEK_RasT_CRaf_CRaf, 2);
            break;
        case 766:
            dJydy[0] = (-1.0*myFOXOnuc + 1.0*yFOXOnuc)/std::pow(sigma_yFOXOnuc, 2);
            break;
        case 767:
            dJydy[0] = (-1.0*myMEKi + 1.0*yMEKi)/std::pow(sigma_yMEKi, 2);
            break;
        case 768:
            dJydy[0] = (-1.0*myMEKi_ppMEK + 1.0*yMEKi_ppMEK)/std::pow(sigma_yMEKi_ppMEK, 2);
            break;
        case 769:
            dJydy[0] = (-1.0*myAKTi + 1.0*yAKTi)/std::pow(sigma_yAKTi, 2);
            break;
        case 770:
            dJydy[0] = (-1.0*myAKTi_AKT + 1.0*yAKTi_AKT)/std::pow(sigma_yAKTi_AKT, 2);
            break;
        case 771:
            dJydy[0] = (-1.0*mymT + 1.0*ymT)/std::pow(sigma_ymT, 2);
            break;
        case 772:
            dJydy[0] = (-1.0*myEIF4E_mT + 1.0*yEIF4E_mT)/std::pow(sigma_yEIF4E_mT, 2);
            break;
        case 773:
            dJydy[0] = (-1.0*mym_TP53 + 1.0*ym_TP53)/std::pow(sigma_ym_TP53, 2);
            break;
        case 774:
            dJydy[0] = (-1.0*mym_MDM2 + 1.0*ym_MDM2)/std::pow(sigma_ym_MDM2, 2);
            break;
        case 775:
            dJydy[0] = (-1.0*mym_PPM1D + 1.0*ym_PPM1D)/std::pow(sigma_ym_PPM1D, 2);
            break;
        case 776:
            dJydy[0] = (-1.0*mym_ATM + 1.0*ym_ATM)/std::pow(sigma_ym_ATM, 2);
            break;
        case 777:
            dJydy[0] = (-1.0*mym_ATR + 1.0*ym_ATR)/std::pow(sigma_ym_ATR, 2);
            break;
        case 778:
            dJydy[0] = (-1.0*mym_RB1 + 1.0*ym_RB1)/std::pow(sigma_ym_RB1, 2);
            break;
        case 779:
            dJydy[0] = (-1.0*mym_E2F1 + 1.0*ym_E2F1)/std::pow(sigma_ym_E2F1, 2);
            break;
        case 780:
            dJydy[0] = (-1.0*mym_E2F2 + 1.0*ym_E2F2)/std::pow(sigma_ym_E2F2, 2);
            break;
        case 781:
            dJydy[0] = (-1.0*mym_E2F3 + 1.0*ym_E2F3)/std::pow(sigma_ym_E2F3, 2);
            break;
        case 782:
            dJydy[0] = (-1.0*mym_CCND1 + 1.0*ym_CCND1)/std::pow(sigma_ym_CCND1, 2);
            break;
        case 783:
            dJydy[0] = (-1.0*mym_CCND2 + 1.0*ym_CCND2)/std::pow(sigma_ym_CCND2, 2);
            break;
        case 784:
            dJydy[0] = (-1.0*mym_CCND3 + 1.0*ym_CCND3)/std::pow(sigma_ym_CCND3, 2);
            break;
        case 785:
            dJydy[0] = (-1.0*mym_CCNE1 + 1.0*ym_CCNE1)/std::pow(sigma_ym_CCNE1, 2);
            break;
        case 786:
            dJydy[0] = (-1.0*mym_CCNE2 + 1.0*ym_CCNE2)/std::pow(sigma_ym_CCNE2, 2);
            break;
        case 787:
            dJydy[0] = (-1.0*mym_SKP2 + 1.0*ym_SKP2)/std::pow(sigma_ym_SKP2, 2);
            break;
        case 788:
            dJydy[0] = (-1.0*mym_CDC25A + 1.0*ym_CDC25A)/std::pow(sigma_ym_CDC25A, 2);
            break;
        case 789:
            dJydy[0] = (-1.0*mym_CDC25B + 1.0*ym_CDC25B)/std::pow(sigma_ym_CDC25B, 2);
            break;
        case 790:
            dJydy[0] = (-1.0*mym_CDC25C + 1.0*ym_CDC25C)/std::pow(sigma_ym_CDC25C, 2);
            break;
        case 791:
            dJydy[0] = (-1.0*mym_CCNA2 + 1.0*ym_CCNA2)/std::pow(sigma_ym_CCNA2, 2);
            break;
        case 792:
            dJydy[0] = (-1.0*mym_CDKN1B + 1.0*ym_CDKN1B)/std::pow(sigma_ym_CDKN1B, 2);
            break;
        case 793:
            dJydy[0] = (-1.0*mym_CDH1 + 1.0*ym_CDH1)/std::pow(sigma_ym_CDH1, 2);
            break;
        case 794:
            dJydy[0] = (-1.0*mym_CCNB1 + 1.0*ym_CCNB1)/std::pow(sigma_ym_CCNB1, 2);
            break;
        case 795:
            dJydy[0] = (-1.0*mym_CDC20 + 1.0*ym_CDC20)/std::pow(sigma_ym_CDC20, 2);
            break;
        case 796:
            dJydy[0] = (-1.0*mym_WEE1 + 1.0*ym_WEE1)/std::pow(sigma_ym_WEE1, 2);
            break;
        case 797:
            dJydy[0] = (-1.0*mym_CHEK1 + 1.0*ym_CHEK1)/std::pow(sigma_ym_CHEK1, 2);
            break;
        case 798:
            dJydy[0] = (-1.0*mym_CDKN1A + 1.0*ym_CDKN1A)/std::pow(sigma_ym_CDKN1A, 2);
            break;
        case 799:
            dJydy[0] = (-1.0*mym_CDK1 + 1.0*ym_CDK1)/std::pow(sigma_ym_CDK1, 2);
            break;
        case 800:
            dJydy[0] = (-1.0*mym_CDK2 + 1.0*ym_CDK2)/std::pow(sigma_ym_CDK2, 2);
            break;
        case 801:
            dJydy[0] = (-1.0*mym_CDK4 + 1.0*ym_CDK4)/std::pow(sigma_ym_CDK4, 2);
            break;
        case 802:
            dJydy[0] = (-1.0*mym_CDK6 + 1.0*ym_CDK6)/std::pow(sigma_ym_CDK6, 2);
            break;
        case 803:
            dJydy[0] = (-1.0*mym_TNFSF10 + 1.0*ym_TNFSF10)/std::pow(sigma_ym_TNFSF10, 2);
            break;
        case 804:
            dJydy[0] = (-1.0*mym_TNFRSF10A + 1.0*ym_TNFRSF10A)/std::pow(sigma_ym_TNFRSF10A, 2);
            break;
        case 805:
            dJydy[0] = (-1.0*mym_TNFRSF10B + 1.0*ym_TNFRSF10B)/std::pow(sigma_ym_TNFRSF10B, 2);
            break;
        case 806:
            dJydy[0] = (-1.0*mym_CFLAR + 1.0*ym_CFLAR)/std::pow(sigma_ym_CFLAR, 2);
            break;
        case 807:
            dJydy[0] = (-1.0*mym_CASP8 + 1.0*ym_CASP8)/std::pow(sigma_ym_CASP8, 2);
            break;
        case 808:
            dJydy[0] = (-1.0*mym_CASP10 + 1.0*ym_CASP10)/std::pow(sigma_ym_CASP10, 2);
            break;
        case 809:
            dJydy[0] = (-1.0*mym_BFAR + 1.0*ym_BFAR)/std::pow(sigma_ym_BFAR, 2);
            break;
        case 810:
            dJydy[0] = (-1.0*mym_CASP3 + 1.0*ym_CASP3)/std::pow(sigma_ym_CASP3, 2);
            break;
        case 811:
            dJydy[0] = (-1.0*mym_CASP7 + 1.0*ym_CASP7)/std::pow(sigma_ym_CASP7, 2);
            break;
        case 812:
            dJydy[0] = (-1.0*mym_CASP6 + 1.0*ym_CASP6)/std::pow(sigma_ym_CASP6, 2);
            break;
        case 813:
            dJydy[0] = (-1.0*mym_XIAP + 1.0*ym_XIAP)/std::pow(sigma_ym_XIAP, 2);
            break;
        case 814:
            dJydy[0] = (-1.0*mym_PARP1 + 1.0*ym_PARP1)/std::pow(sigma_ym_PARP1, 2);
            break;
        case 815:
            dJydy[0] = (-1.0*mym_BID + 1.0*ym_BID)/std::pow(sigma_ym_BID, 2);
            break;
        case 816:
            dJydy[0] = (-1.0*mym_BCL2 + 1.0*ym_BCL2)/std::pow(sigma_ym_BCL2, 2);
            break;
        case 817:
            dJydy[0] = (-1.0*mym_BCL2L1 + 1.0*ym_BCL2L1)/std::pow(sigma_ym_BCL2L1, 2);
            break;
        case 818:
            dJydy[0] = (-1.0*mym_MCL1 + 1.0*ym_MCL1)/std::pow(sigma_ym_MCL1, 2);
            break;
        case 819:
            dJydy[0] = (-1.0*mym_BAX + 1.0*ym_BAX)/std::pow(sigma_ym_BAX, 2);
            break;
        case 820:
            dJydy[0] = (-1.0*mym_CYCS + 1.0*ym_CYCS)/std::pow(sigma_ym_CYCS, 2);
            break;
        case 821:
            dJydy[0] = (-1.0*mym_DIABLO + 1.0*ym_DIABLO)/std::pow(sigma_ym_DIABLO, 2);
            break;
        case 822:
            dJydy[0] = (-1.0*mym_APAF1 + 1.0*ym_APAF1)/std::pow(sigma_ym_APAF1, 2);
            break;
        case 823:
            dJydy[0] = (-1.0*mym_CASP9 + 1.0*ym_CASP9)/std::pow(sigma_ym_CASP9, 2);
            break;
        case 824:
            dJydy[0] = (-1.0*mym_BAD + 1.0*ym_BAD)/std::pow(sigma_ym_BAD, 2);
            break;
        case 825:
            dJydy[0] = (-1.0*mym_BBC3 + 1.0*ym_BBC3)/std::pow(sigma_ym_BBC3, 2);
            break;
        case 826:
            dJydy[0] = (-1.0*mym_PMAIP1 + 1.0*ym_PMAIP1)/std::pow(sigma_ym_PMAIP1, 2);
            break;
        case 827:
            dJydy[0] = (-1.0*mym_BCL2L11 + 1.0*ym_BCL2L11)/std::pow(sigma_ym_BCL2L11, 2);
            break;
        case 828:
            dJydy[0] = (-1.0*mym_EGF + 1.0*ym_EGF)/std::pow(sigma_ym_EGF, 2);
            break;
        case 829:
            dJydy[0] = (-1.0*mym_NRG1 + 1.0*ym_NRG1)/std::pow(sigma_ym_NRG1, 2);
            break;
        case 830:
            dJydy[0] = (-1.0*mym_EGFR + 1.0*ym_EGFR)/std::pow(sigma_ym_EGFR, 2);
            break;
        case 831:
            dJydy[0] = (-1.0*mym_ERBB2 + 1.0*ym_ERBB2)/std::pow(sigma_ym_ERBB2, 2);
            break;
        case 832:
            dJydy[0] = (-1.0*mym_ERBB3 + 1.0*ym_ERBB3)/std::pow(sigma_ym_ERBB3, 2);
            break;
        case 833:
            dJydy[0] = (-1.0*mym_ERBB4 + 1.0*ym_ERBB4)/std::pow(sigma_ym_ERBB4, 2);
            break;
        case 834:
            dJydy[0] = (-1.0*mym_EGFRvIII + 1.0*ym_EGFRvIII)/std::pow(sigma_ym_EGFRvIII, 2);
            break;
        case 835:
            dJydy[0] = (-1.0*mym_MET + 1.0*ym_MET)/std::pow(sigma_ym_MET, 2);
            break;
        case 836:
            dJydy[0] = (-1.0*mym_HGF + 1.0*ym_HGF)/std::pow(sigma_ym_HGF, 2);
            break;
        case 837:
            dJydy[0] = (-1.0*mym_PDGFRA + 1.0*ym_PDGFRA)/std::pow(sigma_ym_PDGFRA, 2);
            break;
        case 838:
            dJydy[0] = (-1.0*mym_PDGFRB + 1.0*ym_PDGFRB)/std::pow(sigma_ym_PDGFRB, 2);
            break;
        case 839:
            dJydy[0] = (-1.0*mym_PDGFB + 1.0*ym_PDGFB)/std::pow(sigma_ym_PDGFB, 2);
            break;
        case 840:
            dJydy[0] = (-1.0*mym_SPRY2 + 1.0*ym_SPRY2)/std::pow(sigma_ym_SPRY2, 2);
            break;
        case 841:
            dJydy[0] = (-1.0*mym_CBL + 1.0*ym_CBL)/std::pow(sigma_ym_CBL, 2);
            break;
        case 842:
            dJydy[0] = (-1.0*mym_GRB2 + 1.0*ym_GRB2)/std::pow(sigma_ym_GRB2, 2);
            break;
        case 843:
            dJydy[0] = (-1.0*mym_PLCG1 + 1.0*ym_PLCG1)/std::pow(sigma_ym_PLCG1, 2);
            break;
        case 844:
            dJydy[0] = (-1.0*mym_PLCG2 + 1.0*ym_PLCG2)/std::pow(sigma_ym_PLCG2, 2);
            break;
        case 845:
            dJydy[0] = (-1.0*mym_PIK3CA + 1.0*ym_PIK3CA)/std::pow(sigma_ym_PIK3CA, 2);
            break;
        case 846:
            dJydy[0] = (-1.0*mym_PIK3CB + 1.0*ym_PIK3CB)/std::pow(sigma_ym_PIK3CB, 2);
            break;
        case 847:
            dJydy[0] = (-1.0*mym_PIK3CG + 1.0*ym_PIK3CG)/std::pow(sigma_ym_PIK3CG, 2);
            break;
        case 848:
            dJydy[0] = (-1.0*mym_PIK3CD + 1.0*ym_PIK3CD)/std::pow(sigma_ym_PIK3CD, 2);
            break;
        case 849:
            dJydy[0] = (-1.0*mym_PIK3R1 + 1.0*ym_PIK3R1)/std::pow(sigma_ym_PIK3R1, 2);
            break;
        case 850:
            dJydy[0] = (-1.0*mym_PIK3R2 + 1.0*ym_PIK3R2)/std::pow(sigma_ym_PIK3R2, 2);
            break;
        case 851:
            dJydy[0] = (-1.0*mym_PIK3R3 + 1.0*ym_PIK3R3)/std::pow(sigma_ym_PIK3R3, 2);
            break;
        case 852:
            dJydy[0] = (-1.0*mym_PIK3R4 + 1.0*ym_PIK3R4)/std::pow(sigma_ym_PIK3R4, 2);
            break;
        case 853:
            dJydy[0] = (-1.0*mym_PIK3C2A + 1.0*ym_PIK3C2A)/std::pow(sigma_ym_PIK3C2A, 2);
            break;
        case 854:
            dJydy[0] = (-1.0*mym_RASGRP1 + 1.0*ym_RASGRP1)/std::pow(sigma_ym_RASGRP1, 2);
            break;
        case 855:
            dJydy[0] = (-1.0*mym_RASGRP3 + 1.0*ym_RASGRP3)/std::pow(sigma_ym_RASGRP3, 2);
            break;
        case 856:
            dJydy[0] = (-1.0*mym_NRAS + 1.0*ym_NRAS)/std::pow(sigma_ym_NRAS, 2);
            break;
        case 857:
            dJydy[0] = (-1.0*mym_KRAS + 1.0*ym_KRAS)/std::pow(sigma_ym_KRAS, 2);
            break;
        case 858:
            dJydy[0] = (-1.0*mym_HRAS + 1.0*ym_HRAS)/std::pow(sigma_ym_HRAS, 2);
            break;
        case 859:
            dJydy[0] = (-1.0*mym_NF1 + 1.0*ym_NF1)/std::pow(sigma_ym_NF1, 2);
            break;
        case 860:
            dJydy[0] = (-1.0*mym_RAF1 + 1.0*ym_RAF1)/std::pow(sigma_ym_RAF1, 2);
            break;
        case 861:
            dJydy[0] = (-1.0*mym_BRAF + 1.0*ym_BRAF)/std::pow(sigma_ym_BRAF, 2);
            break;
        case 862:
            dJydy[0] = (-1.0*mym_MAP2K1 + 1.0*ym_MAP2K1)/std::pow(sigma_ym_MAP2K1, 2);
            break;
        case 863:
            dJydy[0] = (-1.0*mym_MAP2K2 + 1.0*ym_MAP2K2)/std::pow(sigma_ym_MAP2K2, 2);
            break;
        case 864:
            dJydy[0] = (-1.0*mym_DUSP6 + 1.0*ym_DUSP6)/std::pow(sigma_ym_DUSP6, 2);
            break;
        case 865:
            dJydy[0] = (-1.0*mym_RPS6KA1 + 1.0*ym_RPS6KA1)/std::pow(sigma_ym_RPS6KA1, 2);
            break;
        case 866:
            dJydy[0] = (-1.0*mym_RPS6KA2 + 1.0*ym_RPS6KA2)/std::pow(sigma_ym_RPS6KA2, 2);
            break;
        case 867:
            dJydy[0] = (-1.0*mym_RPS6KA3 + 1.0*ym_RPS6KA3)/std::pow(sigma_ym_RPS6KA3, 2);
            break;
        case 868:
            dJydy[0] = (-1.0*mym_RPS6KA4 + 1.0*ym_RPS6KA4)/std::pow(sigma_ym_RPS6KA4, 2);
            break;
        case 869:
            dJydy[0] = (-1.0*mym_DUSP1 + 1.0*ym_DUSP1)/std::pow(sigma_ym_DUSP1, 2);
            break;
        case 870:
            dJydy[0] = (-1.0*mym_FOS + 1.0*ym_FOS)/std::pow(sigma_ym_FOS, 2);
            break;
        case 871:
            dJydy[0] = (-1.0*mym_JUN + 1.0*ym_JUN)/std::pow(sigma_ym_JUN, 2);
            break;
        case 872:
            dJydy[0] = (-1.0*mym_MYC + 1.0*ym_MYC)/std::pow(sigma_ym_MYC, 2);
            break;
        case 873:
            dJydy[0] = (-1.0*mym_CTNNB1 + 1.0*ym_CTNNB1)/std::pow(sigma_ym_CTNNB1, 2);
            break;
        case 874:
            dJydy[0] = (-1.0*mym_PTEN + 1.0*ym_PTEN)/std::pow(sigma_ym_PTEN, 2);
            break;
        case 875:
            dJydy[0] = (-1.0*mym_AKT1 + 1.0*ym_AKT1)/std::pow(sigma_ym_AKT1, 2);
            break;
        case 876:
            dJydy[0] = (-1.0*mym_AKT2 + 1.0*ym_AKT2)/std::pow(sigma_ym_AKT2, 2);
            break;
        case 877:
            dJydy[0] = (-1.0*mym_PDPK1 + 1.0*ym_PDPK1)/std::pow(sigma_ym_PDPK1, 2);
            break;
        case 878:
            dJydy[0] = (-1.0*mym_RICTOR + 1.0*ym_RICTOR)/std::pow(sigma_ym_RICTOR, 2);
            break;
        case 879:
            dJydy[0] = (-1.0*mym_MTOR + 1.0*ym_MTOR)/std::pow(sigma_ym_MTOR, 2);
            break;
        case 880:
            dJydy[0] = (-1.0*mym_GSK3B + 1.0*ym_GSK3B)/std::pow(sigma_ym_GSK3B, 2);
            break;
        case 881:
            dJydy[0] = (-1.0*mym_TSC1 + 1.0*ym_TSC1)/std::pow(sigma_ym_TSC1, 2);
            break;
        case 882:
            dJydy[0] = (-1.0*mym_TSC2 + 1.0*ym_TSC2)/std::pow(sigma_ym_TSC2, 2);
            break;
        case 883:
            dJydy[0] = (-1.0*mym_PRKCA + 1.0*ym_PRKCA)/std::pow(sigma_ym_PRKCA, 2);
            break;
        case 884:
            dJydy[0] = (-1.0*mym_PRKCB + 1.0*ym_PRKCB)/std::pow(sigma_ym_PRKCB, 2);
            break;
        case 885:
            dJydy[0] = (-1.0*mym_PRKCG + 1.0*ym_PRKCG)/std::pow(sigma_ym_PRKCG, 2);
            break;
        case 886:
            dJydy[0] = (-1.0*mym_PRKCD + 1.0*ym_PRKCD)/std::pow(sigma_ym_PRKCD, 2);
            break;
        case 887:
            dJydy[0] = (-1.0*mym_PEBP1 + 1.0*ym_PEBP1)/std::pow(sigma_ym_PEBP1, 2);
            break;
        case 888:
            dJydy[0] = (-1.0*mym_MAPK1 + 1.0*ym_MAPK1)/std::pow(sigma_ym_MAPK1, 2);
            break;
        case 889:
            dJydy[0] = (-1.0*mym_MAPK3 + 1.0*ym_MAPK3)/std::pow(sigma_ym_MAPK3, 2);
            break;
        case 890:
            dJydy[0] = (-1.0*mym_FOXO3 + 1.0*ym_FOXO3)/std::pow(sigma_ym_FOXO3, 2);
            break;
        case 891:
            dJydy[0] = (-1.0*mym_RHEB + 1.0*ym_RHEB)/std::pow(sigma_ym_RHEB, 2);
            break;
        case 892:
            dJydy[0] = (-1.0*mym_RPTOR + 1.0*ym_RPTOR)/std::pow(sigma_ym_RPTOR, 2);
            break;
        case 893:
            dJydy[0] = (-1.0*mym_RPS6KB1 + 1.0*ym_RPS6KB1)/std::pow(sigma_ym_RPS6KB1, 2);
            break;
        case 894:
            dJydy[0] = (-1.0*mym_RPS6KB2 + 1.0*ym_RPS6KB2)/std::pow(sigma_ym_RPS6KB2, 2);
            break;
        case 895:
            dJydy[0] = (-1.0*mym_EIF4EBP1 + 1.0*ym_EIF4EBP1)/std::pow(sigma_ym_EIF4EBP1, 2);
            break;
        case 896:
            dJydy[0] = (-1.0*mym_SOS1 + 1.0*ym_SOS1)/std::pow(sigma_ym_SOS1, 2);
            break;
        case 897:
            dJydy[0] = (-1.0*mym_CDKN2A + 1.0*ym_CDKN2A)/std::pow(sigma_ym_CDKN2A, 2);
            break;
        case 898:
            dJydy[0] = (-1.0*mym_MDM4 + 1.0*ym_MDM4)/std::pow(sigma_ym_MDM4, 2);
            break;
        case 899:
            dJydy[0] = (-1.0*mym_FGFR1 + 1.0*ym_FGFR1)/std::pow(sigma_ym_FGFR1, 2);
            break;
        case 900:
            dJydy[0] = (-1.0*mym_FGFR2 + 1.0*ym_FGFR2)/std::pow(sigma_ym_FGFR2, 2);
            break;
        case 901:
            dJydy[0] = (-1.0*mym_FGF1 + 1.0*ym_FGF1)/std::pow(sigma_ym_FGF1, 2);
            break;
        case 902:
            dJydy[0] = (-1.0*mym_FGF2 + 1.0*ym_FGF2)/std::pow(sigma_ym_FGF2, 2);
            break;
        case 903:
            dJydy[0] = (-1.0*mym_EIF4E + 1.0*ym_EIF4E)/std::pow(sigma_ym_EIF4E, 2);
            break;
        case 904:
            dJydy[0] = (-1.0*mym_IRS1 + 1.0*ym_IRS1)/std::pow(sigma_ym_IRS1, 2);
            break;
        case 905:
            dJydy[0] = (-1.0*mym_IRS2 + 1.0*ym_IRS2)/std::pow(sigma_ym_IRS2, 2);
            break;
        case 906:
            dJydy[0] = (-1.0*mym_IGF1 + 1.0*ym_IGF1)/std::pow(sigma_ym_IGF1, 2);
            break;
        case 907:
            dJydy[0] = (-1.0*mym_IGF2 + 1.0*ym_IGF2)/std::pow(sigma_ym_IGF2, 2);
            break;
        case 908:
            dJydy[0] = (-1.0*mym_IGF1R + 1.0*ym_IGF1R)/std::pow(sigma_ym_IGF1R, 2);
            break;
        case 909:
            dJydy[0] = (-1.0*mym_MSH6 + 1.0*ym_MSH6)/std::pow(sigma_ym_MSH6, 2);
            break;
        case 910:
            dJydy[0] = (-1.0*mym_BRCA2 + 1.0*ym_BRCA2)/std::pow(sigma_ym_BRCA2, 2);
            break;
        case 911:
            dJydy[0] = (-1.0*mym_MGMT + 1.0*ym_MGMT)/std::pow(sigma_ym_MGMT, 2);
            break;
        case 912:
            dJydy[0] = (-1.0*mym_INSR + 1.0*ym_INSR)/std::pow(sigma_ym_INSR, 2);
            break;
        case 913:
            dJydy[0] = (-1.0*mym_INS + 1.0*ym_INS)/std::pow(sigma_ym_INS, 2);
            break;
        case 914:
            dJydy[0] = (-1.0*myCytoplasm + 1.0*yCytoplasm)/std::pow(sigma_yCytoplasm, 2);
            break;
        case 915:
            dJydy[0] = (-1.0*myExtracellular + 1.0*yExtracellular)/std::pow(sigma_yExtracellular, 2);
            break;
        case 916:
            dJydy[0] = (-1.0*myNucleus + 1.0*yNucleus)/std::pow(sigma_yNucleus, 2);
            break;
        case 917:
            dJydy[0] = (-1.0*myMitochondrion + 1.0*yMitochondrion)/std::pow(sigma_yMitochondrion, 2);
            break;
    }
}

} // namespace amici
} // namespace model_SPARCED_tutorial