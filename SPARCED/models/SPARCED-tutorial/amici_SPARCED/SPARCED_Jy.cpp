#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <array>

#include "y.h"
#include "my.h"
#include "p.h"
#include "k.h"
#include "sigmay.h"

namespace amici {
namespace model_SPARCED {

void Jy_SPARCED(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yRibosome, 2)) + 0.5*std::pow(-myRibosome + yRibosome, 2)/std::pow(sigma_yRibosome, 2);
            break;
        case 1:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yp53inac, 2)) + 0.5*std::pow(-myp53inac + yp53inac, 2)/std::pow(sigma_yp53inac, 2);
            break;
        case 2:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yp53ac, 2)) + 0.5*std::pow(-myp53ac + yp53ac, 2)/std::pow(sigma_yp53ac, 2);
            break;
        case 3:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMDM2, 2)) + 0.5*std::pow(-myMDM2 + yMDM2, 2)/std::pow(sigma_yMDM2, 2);
            break;
        case 4:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yWip1, 2)) + 0.5*std::pow(-myWip1 + yWip1, 2)/std::pow(sigma_yWip1, 2);
            break;
        case 5:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yATMP, 2)) + 0.5*std::pow(-myATMP + yATMP, 2)/std::pow(sigma_yATMP, 2);
            break;
        case 6:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yATRac, 2)) + 0.5*std::pow(-myATRac + yATRac, 2)/std::pow(sigma_yATRac, 2);
            break;
        case 7:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMDM2product1, 2)) + 0.5*std::pow(-myMDM2product1 + yMDM2product1, 2)/std::pow(sigma_yMDM2product1, 2);
            break;
        case 8:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMDM2product2, 2)) + 0.5*std::pow(-myMDM2product2 + yMDM2product2, 2)/std::pow(sigma_yMDM2product2, 2);
            break;
        case 9:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMDM2product3, 2)) + 0.5*std::pow(-myMDM2product3 + yMDM2product3, 2)/std::pow(sigma_yMDM2product3, 2);
            break;
        case 10:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMDM2product4, 2)) + 0.5*std::pow(-myMDM2product4 + yMDM2product4, 2)/std::pow(sigma_yMDM2product4, 2);
            break;
        case 11:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMDM2product5, 2)) + 0.5*std::pow(-myMDM2product5 + yMDM2product5, 2)/std::pow(sigma_yMDM2product5, 2);
            break;
        case 12:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMDM2product6, 2)) + 0.5*std::pow(-myMDM2product6 + yMDM2product6, 2)/std::pow(sigma_yMDM2product6, 2);
            break;
        case 13:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMDM2product7, 2)) + 0.5*std::pow(-myMDM2product7 + yMDM2product7, 2)/std::pow(sigma_yMDM2product7, 2);
            break;
        case 14:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMDM2product8, 2)) + 0.5*std::pow(-myMDM2product8 + yMDM2product8, 2)/std::pow(sigma_yMDM2product8, 2);
            break;
        case 15:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMDM2product9, 2)) + 0.5*std::pow(-myMDM2product9 + yMDM2product9, 2)/std::pow(sigma_yMDM2product9, 2);
            break;
        case 16:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMDM2pro, 2)) + 0.5*std::pow(-myMDM2pro + yMDM2pro, 2)/std::pow(sigma_yMDM2pro, 2);
            break;
        case 17:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yWip1product1, 2)) + 0.5*std::pow(-myWip1product1 + yWip1product1, 2)/std::pow(sigma_yWip1product1, 2);
            break;
        case 18:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yWip1product2, 2)) + 0.5*std::pow(-myWip1product2 + yWip1product2, 2)/std::pow(sigma_yWip1product2, 2);
            break;
        case 19:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yWip1product3, 2)) + 0.5*std::pow(-myWip1product3 + yWip1product3, 2)/std::pow(sigma_yWip1product3, 2);
            break;
        case 20:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yWip1product4, 2)) + 0.5*std::pow(-myWip1product4 + yWip1product4, 2)/std::pow(sigma_yWip1product4, 2);
            break;
        case 21:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yWip1product5, 2)) + 0.5*std::pow(-myWip1product5 + yWip1product5, 2)/std::pow(sigma_yWip1product5, 2);
            break;
        case 22:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yWip1product6, 2)) + 0.5*std::pow(-myWip1product6 + yWip1product6, 2)/std::pow(sigma_yWip1product6, 2);
            break;
        case 23:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yWip1product7, 2)) + 0.5*std::pow(-myWip1product7 + yWip1product7, 2)/std::pow(sigma_yWip1product7, 2);
            break;
        case 24:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yWip1product8, 2)) + 0.5*std::pow(-myWip1product8 + yWip1product8, 2)/std::pow(sigma_yWip1product8, 2);
            break;
        case 25:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yWip1product9, 2)) + 0.5*std::pow(-myWip1product9 + yWip1product9, 2)/std::pow(sigma_yWip1product9, 2);
            break;
        case 26:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yWip1pro, 2)) + 0.5*std::pow(-myWip1pro + yWip1pro, 2)/std::pow(sigma_yWip1pro, 2);
            break;
        case 27:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yBRCA2, 2)) + 0.5*std::pow(-myBRCA2 + yBRCA2, 2)/std::pow(sigma_yBRCA2, 2);
            break;
        case 28:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMSH6, 2)) + 0.5*std::pow(-myMSH6 + yMSH6, 2)/std::pow(sigma_yMSH6, 2);
            break;
        case 29:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMGMT, 2)) + 0.5*std::pow(-myMGMT + yMGMT, 2)/std::pow(sigma_yMGMT, 2);
            break;
        case 30:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ydamageDSB, 2)) + 0.5*std::pow(-mydamageDSB + ydamageDSB, 2)/std::pow(sigma_ydamageDSB, 2);
            break;
        case 31:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ydamageSSB, 2)) + 0.5*std::pow(-mydamageSSB + ydamageSSB, 2)/std::pow(sigma_ydamageSSB, 2);
            break;
        case 32:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yppAKT_MDM2, 2)) + 0.5*std::pow(-myppAKT_MDM2 + yppAKT_MDM2, 2)/std::pow(sigma_yppAKT_MDM2, 2);
            break;
        case 33:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypMDM2, 2)) + 0.5*std::pow(-mypMDM2 + ypMDM2, 2)/std::pow(sigma_ypMDM2, 2);
            break;
        case 34:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yARF, 2)) + 0.5*std::pow(-myARF + yARF, 2)/std::pow(sigma_yARF, 2);
            break;
        case 35:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMDM4, 2)) + 0.5*std::pow(-myMDM4 + yMDM4, 2)/std::pow(sigma_yMDM4, 2);
            break;
        case 36:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yp53ac_MDM4, 2)) + 0.5*std::pow(-myp53ac_MDM4 + yp53ac_MDM4, 2)/std::pow(sigma_yp53ac_MDM4, 2);
            break;
        case 37:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yATMinac, 2)) + 0.5*std::pow(-myATMinac + yATMinac, 2)/std::pow(sigma_yATMinac, 2);
            break;
        case 38:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yATRinac, 2)) + 0.5*std::pow(-myATRinac + yATRinac, 2)/std::pow(sigma_yATRinac, 2);
            break;
        case 39:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypRB, 2)) + 0.5*std::pow(-mypRB + ypRB, 2)/std::pow(sigma_ypRB, 2);
            break;
        case 40:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypRBp, 2)) + 0.5*std::pow(-mypRBp + ypRBp, 2)/std::pow(sigma_ypRBp, 2);
            break;
        case 41:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypRBpp, 2)) + 0.5*std::pow(-mypRBpp + ypRBpp, 2)/std::pow(sigma_ypRBpp, 2);
            break;
        case 42:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yE2F, 2)) + 0.5*std::pow(-myE2F + yE2F, 2)/std::pow(sigma_yE2F, 2);
            break;
        case 43:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yCd, 2)) + 0.5*std::pow(-myCd + yCd, 2)/std::pow(sigma_yCd, 2);
            break;
        case 44:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMdi, 2)) + 0.5*std::pow(-myMdi + yMdi, 2)/std::pow(sigma_yMdi, 2);
            break;
        case 45:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMd, 2)) + 0.5*std::pow(-myMd + yMd, 2)/std::pow(sigma_yMd, 2);
            break;
        case 46:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMdp27, 2)) + 0.5*std::pow(-myMdp27 + yMdp27, 2)/std::pow(sigma_yMdp27, 2);
            break;
        case 47:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yCe, 2)) + 0.5*std::pow(-myCe + yCe, 2)/std::pow(sigma_yCe, 2);
            break;
        case 48:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMei, 2)) + 0.5*std::pow(-myMei + yMei, 2)/std::pow(sigma_yMei, 2);
            break;
        case 49:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMe, 2)) + 0.5*std::pow(-myMe + yMe, 2)/std::pow(sigma_yMe, 2);
            break;
        case 50:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ySkp2, 2)) + 0.5*std::pow(-mySkp2 + ySkp2, 2)/std::pow(sigma_ySkp2, 2);
            break;
        case 51:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMep27, 2)) + 0.5*std::pow(-myMep27 + yMep27, 2)/std::pow(sigma_yMep27, 2);
            break;
        case 52:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yPe, 2)) + 0.5*std::pow(-myPe + yPe, 2)/std::pow(sigma_yPe, 2);
            break;
        case 53:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yPai, 2)) + 0.5*std::pow(-myPai + yPai, 2)/std::pow(sigma_yPai, 2);
            break;
        case 54:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yPei, 2)) + 0.5*std::pow(-myPei + yPei, 2)/std::pow(sigma_yPei, 2);
            break;
        case 55:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yPbi, 2)) + 0.5*std::pow(-myPbi + yPbi, 2)/std::pow(sigma_yPbi, 2);
            break;
        case 56:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yCa, 2)) + 0.5*std::pow(-myCa + yCa, 2)/std::pow(sigma_yCa, 2);
            break;
        case 57:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMai, 2)) + 0.5*std::pow(-myMai + yMai, 2)/std::pow(sigma_yMai, 2);
            break;
        case 58:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMa, 2)) + 0.5*std::pow(-myMa + yMa, 2)/std::pow(sigma_yMa, 2);
            break;
        case 59:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMap27, 2)) + 0.5*std::pow(-myMap27 + yMap27, 2)/std::pow(sigma_yMap27, 2);
            break;
        case 60:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yp27, 2)) + 0.5*std::pow(-myp27 + yp27, 2)/std::pow(sigma_yp27, 2);
            break;
        case 61:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yCdh1i, 2)) + 0.5*std::pow(-myCdh1i + yCdh1i, 2)/std::pow(sigma_yCdh1i, 2);
            break;
        case 62:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yCdh1a, 2)) + 0.5*std::pow(-myCdh1a + yCdh1a, 2)/std::pow(sigma_yCdh1a, 2);
            break;
        case 63:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yE2Fp, 2)) + 0.5*std::pow(-myE2Fp + yE2Fp, 2)/std::pow(sigma_yE2Fp, 2);
            break;
        case 64:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yp27p, 2)) + 0.5*std::pow(-myp27p + yp27p, 2)/std::pow(sigma_yp27p, 2);
            break;
        case 65:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yPa, 2)) + 0.5*std::pow(-myPa + yPa, 2)/std::pow(sigma_yPa, 2);
            break;
        case 66:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yCb, 2)) + 0.5*std::pow(-myCb + yCb, 2)/std::pow(sigma_yCb, 2);
            break;
        case 67:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMbi, 2)) + 0.5*std::pow(-myMbi + yMbi, 2)/std::pow(sigma_yMbi, 2);
            break;
        case 68:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMb, 2)) + 0.5*std::pow(-myMb + yMb, 2)/std::pow(sigma_yMb, 2);
            break;
        case 69:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yCdc20i, 2)) + 0.5*std::pow(-myCdc20i + yCdc20i, 2)/std::pow(sigma_yCdc20i, 2);
            break;
        case 70:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yCdc20a, 2)) + 0.5*std::pow(-myCdc20a + yCdc20a, 2)/std::pow(sigma_yCdc20a, 2);
            break;
        case 71:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yPb, 2)) + 0.5*std::pow(-myPb + yPb, 2)/std::pow(sigma_yPb, 2);
            break;
        case 72:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yWee1, 2)) + 0.5*std::pow(-myWee1 + yWee1, 2)/std::pow(sigma_yWee1, 2);
            break;
        case 73:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yWee1p, 2)) + 0.5*std::pow(-myWee1p + yWee1p, 2)/std::pow(sigma_yWee1p, 2);
            break;
        case 74:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMbp27, 2)) + 0.5*std::pow(-myMbp27 + yMbp27, 2)/std::pow(sigma_yMbp27, 2);
            break;
        case 75:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yChk1, 2)) + 0.5*std::pow(-myChk1 + yChk1, 2)/std::pow(sigma_yChk1, 2);
            break;
        case 76:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypRBc1, 2)) + 0.5*std::pow(-mypRBc1 + ypRBc1, 2)/std::pow(sigma_ypRBc1, 2);
            break;
        case 77:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypRBc2, 2)) + 0.5*std::pow(-mypRBc2 + ypRBc2, 2)/std::pow(sigma_ypRBc2, 2);
            break;
        case 78:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yp21, 2)) + 0.5*std::pow(-myp21 + yp21, 2)/std::pow(sigma_yp21, 2);
            break;
        case 79:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMdp21, 2)) + 0.5*std::pow(-myMdp21 + yMdp21, 2)/std::pow(sigma_yMdp21, 2);
            break;
        case 80:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMep21, 2)) + 0.5*std::pow(-myMep21 + yMep21, 2)/std::pow(sigma_yMep21, 2);
            break;
        case 81:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMap21, 2)) + 0.5*std::pow(-myMap21 + yMap21, 2)/std::pow(sigma_yMap21, 2);
            break;
        case 82:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMbp21, 2)) + 0.5*std::pow(-myMbp21 + yMbp21, 2)/std::pow(sigma_yMbp21, 2);
            break;
        case 83:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yL, 2)) + 0.5*std::pow(-myL + yL, 2)/std::pow(sigma_yL, 2);
            break;
        case 84:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yR, 2)) + 0.5*std::pow(-myR + yR, 2)/std::pow(sigma_yR, 2);
            break;
        case 85:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yL_R, 2)) + 0.5*std::pow(-myL_R + yL_R, 2)/std::pow(sigma_yL_R, 2);
            break;
        case 86:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yRactive, 2)) + 0.5*std::pow(-myRactive + yRactive, 2)/std::pow(sigma_yRactive, 2);
            break;
        case 87:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yflip, 2)) + 0.5*std::pow(-myflip + yflip, 2)/std::pow(sigma_yflip, 2);
            break;
        case 88:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yRactive_flip, 2)) + 0.5*std::pow(-myRactive_flip + yRactive_flip, 2)/std::pow(sigma_yRactive_flip, 2);
            break;
        case 89:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypC8, 2)) + 0.5*std::pow(-mypC8 + ypC8, 2)/std::pow(sigma_ypC8, 2);
            break;
        case 90:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yRactive_pC8, 2)) + 0.5*std::pow(-myRactive_pC8 + yRactive_pC8, 2)/std::pow(sigma_yRactive_pC8, 2);
            break;
        case 91:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yC8, 2)) + 0.5*std::pow(-myC8 + yC8, 2)/std::pow(sigma_yC8, 2);
            break;
        case 92:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yBar, 2)) + 0.5*std::pow(-myBar + yBar, 2)/std::pow(sigma_yBar, 2);
            break;
        case 93:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yC8_Bar, 2)) + 0.5*std::pow(-myC8_Bar + yC8_Bar, 2)/std::pow(sigma_yC8_Bar, 2);
            break;
        case 94:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypC3, 2)) + 0.5*std::pow(-mypC3 + ypC3, 2)/std::pow(sigma_ypC3, 2);
            break;
        case 95:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yC8_pC3, 2)) + 0.5*std::pow(-myC8_pC3 + yC8_pC3, 2)/std::pow(sigma_yC8_pC3, 2);
            break;
        case 96:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yC3, 2)) + 0.5*std::pow(-myC3 + yC3, 2)/std::pow(sigma_yC3, 2);
            break;
        case 97:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypC6, 2)) + 0.5*std::pow(-mypC6 + ypC6, 2)/std::pow(sigma_ypC6, 2);
            break;
        case 98:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yC3_pC6, 2)) + 0.5*std::pow(-myC3_pC6 + yC3_pC6, 2)/std::pow(sigma_yC3_pC6, 2);
            break;
        case 99:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yC6, 2)) + 0.5*std::pow(-myC6 + yC6, 2)/std::pow(sigma_yC6, 2);
            break;
        case 100:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yC6_C8, 2)) + 0.5*std::pow(-myC6_C8 + yC6_C8, 2)/std::pow(sigma_yC6_C8, 2);
            break;
        case 101:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yXIAP, 2)) + 0.5*std::pow(-myXIAP + yXIAP, 2)/std::pow(sigma_yXIAP, 2);
            break;
        case 102:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yC3_XIAP, 2)) + 0.5*std::pow(-myC3_XIAP + yC3_XIAP, 2)/std::pow(sigma_yC3_XIAP, 2);
            break;
        case 103:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yPARP, 2)) + 0.5*std::pow(-myPARP + yPARP, 2)/std::pow(sigma_yPARP, 2);
            break;
        case 104:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yC3_PARP, 2)) + 0.5*std::pow(-myC3_PARP + yC3_PARP, 2)/std::pow(sigma_yC3_PARP, 2);
            break;
        case 105:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ycPARP, 2)) + 0.5*std::pow(-mycPARP + ycPARP, 2)/std::pow(sigma_ycPARP, 2);
            break;
        case 106:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yBid, 2)) + 0.5*std::pow(-myBid + yBid, 2)/std::pow(sigma_yBid, 2);
            break;
        case 107:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yC8_Bid, 2)) + 0.5*std::pow(-myC8_Bid + yC8_Bid, 2)/std::pow(sigma_yC8_Bid, 2);
            break;
        case 108:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ytBid, 2)) + 0.5*std::pow(-mytBid + ytBid, 2)/std::pow(sigma_ytBid, 2);
            break;
        case 109:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yBcl2c, 2)) + 0.5*std::pow(-myBcl2c + yBcl2c, 2)/std::pow(sigma_yBcl2c, 2);
            break;
        case 110:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ytBid_Bcl2c, 2)) + 0.5*std::pow(-mytBid_Bcl2c + ytBid_Bcl2c, 2)/std::pow(sigma_ytBid_Bcl2c, 2);
            break;
        case 111:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yBax, 2)) + 0.5*std::pow(-myBax + yBax, 2)/std::pow(sigma_yBax, 2);
            break;
        case 112:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ytBid_Bax, 2)) + 0.5*std::pow(-mytBid_Bax + ytBid_Bax, 2)/std::pow(sigma_ytBid_Bax, 2);
            break;
        case 113:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yBaxactive, 2)) + 0.5*std::pow(-myBaxactive + yBaxactive, 2)/std::pow(sigma_yBaxactive, 2);
            break;
        case 114:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yBaxm, 2)) + 0.5*std::pow(-myBaxm + yBaxm, 2)/std::pow(sigma_yBaxm, 2);
            break;
        case 115:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yBcl2, 2)) + 0.5*std::pow(-myBcl2 + yBcl2, 2)/std::pow(sigma_yBcl2, 2);
            break;
        case 116:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yBaxm_Bcl2, 2)) + 0.5*std::pow(-myBaxm_Bcl2 + yBaxm_Bcl2, 2)/std::pow(sigma_yBaxm_Bcl2, 2);
            break;
        case 117:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yBax2, 2)) + 0.5*std::pow(-myBax2 + yBax2, 2)/std::pow(sigma_yBax2, 2);
            break;
        case 118:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yBax2_Bcl2, 2)) + 0.5*std::pow(-myBax2_Bcl2 + yBax2_Bcl2, 2)/std::pow(sigma_yBax2_Bcl2, 2);
            break;
        case 119:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yBax4, 2)) + 0.5*std::pow(-myBax4 + yBax4, 2)/std::pow(sigma_yBax4, 2);
            break;
        case 120:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yBax4_Bcl2, 2)) + 0.5*std::pow(-myBax4_Bcl2 + yBax4_Bcl2, 2)/std::pow(sigma_yBax4_Bcl2, 2);
            break;
        case 121:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yM, 2)) + 0.5*std::pow(-myM + yM, 2)/std::pow(sigma_yM, 2);
            break;
        case 122:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yBax4_M, 2)) + 0.5*std::pow(-myBax4_M + yBax4_M, 2)/std::pow(sigma_yBax4_M, 2);
            break;
        case 123:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMactive, 2)) + 0.5*std::pow(-myMactive + yMactive, 2)/std::pow(sigma_yMactive, 2);
            break;
        case 124:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yCytoCm, 2)) + 0.5*std::pow(-myCytoCm + yCytoCm, 2)/std::pow(sigma_yCytoCm, 2);
            break;
        case 125:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMactive_CytoCm, 2)) + 0.5*std::pow(-myMactive_CytoCm + yMactive_CytoCm, 2)/std::pow(sigma_yMactive_CytoCm, 2);
            break;
        case 126:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yCytoCr, 2)) + 0.5*std::pow(-myCytoCr + yCytoCr, 2)/std::pow(sigma_yCytoCr, 2);
            break;
        case 127:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ySmacm, 2)) + 0.5*std::pow(-mySmacm + ySmacm, 2)/std::pow(sigma_ySmacm, 2);
            break;
        case 128:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMactive_Smacm, 2)) + 0.5*std::pow(-myMactive_Smacm + yMactive_Smacm, 2)/std::pow(sigma_yMactive_Smacm, 2);
            break;
        case 129:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ySmacr, 2)) + 0.5*std::pow(-mySmacr + ySmacr, 2)/std::pow(sigma_ySmacr, 2);
            break;
        case 130:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yCytoC, 2)) + 0.5*std::pow(-myCytoC + yCytoC, 2)/std::pow(sigma_yCytoC, 2);
            break;
        case 131:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yApaf, 2)) + 0.5*std::pow(-myApaf + yApaf, 2)/std::pow(sigma_yApaf, 2);
            break;
        case 132:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yCytoC_Apaf, 2)) + 0.5*std::pow(-myCytoC_Apaf + yCytoC_Apaf, 2)/std::pow(sigma_yCytoC_Apaf, 2);
            break;
        case 133:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yApafactive, 2)) + 0.5*std::pow(-myApafactive + yApafactive, 2)/std::pow(sigma_yApafactive, 2);
            break;
        case 134:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypC9, 2)) + 0.5*std::pow(-mypC9 + ypC9, 2)/std::pow(sigma_ypC9, 2);
            break;
        case 135:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yApop, 2)) + 0.5*std::pow(-myApop + yApop, 2)/std::pow(sigma_yApop, 2);
            break;
        case 136:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yApop_C3, 2)) + 0.5*std::pow(-myApop_C3 + yApop_C3, 2)/std::pow(sigma_yApop_C3, 2);
            break;
        case 137:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ySmac, 2)) + 0.5*std::pow(-mySmac + ySmac, 2)/std::pow(sigma_ySmac, 2);
            break;
        case 138:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yApop_XIAP, 2)) + 0.5*std::pow(-myApop_XIAP + yApop_XIAP, 2)/std::pow(sigma_yApop_XIAP, 2);
            break;
        case 139:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ySmac_XIAP, 2)) + 0.5*std::pow(-mySmac_XIAP + ySmac_XIAP, 2)/std::pow(sigma_ySmac_XIAP, 2);
            break;
        case 140:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yC3_Ub, 2)) + 0.5*std::pow(-myC3_Ub + yC3_Ub, 2)/std::pow(sigma_yC3_Ub, 2);
            break;
        case 141:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yBAD, 2)) + 0.5*std::pow(-myBAD + yBAD, 2)/std::pow(sigma_yBAD, 2);
            break;
        case 142:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yPUMA, 2)) + 0.5*std::pow(-myPUMA + yPUMA, 2)/std::pow(sigma_yPUMA, 2);
            break;
        case 143:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yNOXA, 2)) + 0.5*std::pow(-myNOXA + yNOXA, 2)/std::pow(sigma_yNOXA, 2);
            break;
        case 144:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yBcl2c_BAD, 2)) + 0.5*std::pow(-myBcl2c_BAD + yBcl2c_BAD, 2)/std::pow(sigma_yBcl2c_BAD, 2);
            break;
        case 145:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yBcl2c_PUMA, 2)) + 0.5*std::pow(-myBcl2c_PUMA + yBcl2c_PUMA, 2)/std::pow(sigma_yBcl2c_PUMA, 2);
            break;
        case 146:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yBcl2c_NOXA, 2)) + 0.5*std::pow(-myBcl2c_NOXA + yBcl2c_NOXA, 2)/std::pow(sigma_yBcl2c_NOXA, 2);
            break;
        case 147:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yBIM, 2)) + 0.5*std::pow(-myBIM + yBIM, 2)/std::pow(sigma_yBIM, 2);
            break;
        case 148:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yBIM_Bax, 2)) + 0.5*std::pow(-myBIM_Bax + yBIM_Bax, 2)/std::pow(sigma_yBIM_Bax, 2);
            break;
        case 149:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yBcl2c_BIM, 2)) + 0.5*std::pow(-myBcl2c_BIM + yBcl2c_BIM, 2)/std::pow(sigma_yBcl2c_BIM, 2);
            break;
        case 150:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yppERK_BIM, 2)) + 0.5*std::pow(-myppERK_BIM + yppERK_BIM, 2)/std::pow(sigma_yppERK_BIM, 2);
            break;
        case 151:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypBIM, 2)) + 0.5*std::pow(-mypBIM + ypBIM, 2)/std::pow(sigma_ypBIM, 2);
            break;
        case 152:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yppAKT_BAD, 2)) + 0.5*std::pow(-myppAKT_BAD + yppAKT_BAD, 2)/std::pow(sigma_yppAKT_BAD, 2);
            break;
        case 153:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypBAD, 2)) + 0.5*std::pow(-mypBAD + ypBAD, 2)/std::pow(sigma_ypBAD, 2);
            break;
        case 154:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yppERK_BAD, 2)) + 0.5*std::pow(-myppERK_BAD + yppERK_BAD, 2)/std::pow(sigma_yppERK_BAD, 2);
            break;
        case 155:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yE, 2)) + 0.5*std::pow(-myE + yE, 2)/std::pow(sigma_yE, 2);
            break;
        case 156:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yH, 2)) + 0.5*std::pow(-myH + yH, 2)/std::pow(sigma_yH, 2);
            break;
        case 157:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yHGF, 2)) + 0.5*std::pow(-myHGF + yHGF, 2)/std::pow(sigma_yHGF, 2);
            break;
        case 158:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yP, 2)) + 0.5*std::pow(-myP + yP, 2)/std::pow(sigma_yP, 2);
            break;
        case 159:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yF, 2)) + 0.5*std::pow(-myF + yF, 2)/std::pow(sigma_yF, 2);
            break;
        case 160:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yI, 2)) + 0.5*std::pow(-myI + yI, 2)/std::pow(sigma_yI, 2);
            break;
        case 161:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yINS, 2)) + 0.5*std::pow(-myINS + yINS, 2)/std::pow(sigma_yINS, 2);
            break;
        case 162:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yE1, 2)) + 0.5*std::pow(-myE1 + yE1, 2)/std::pow(sigma_yE1, 2);
            break;
        case 163:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE1, 2)) + 0.5*std::pow(-mypE1 + ypE1, 2)/std::pow(sigma_ypE1, 2);
            break;
        case 164:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yE2, 2)) + 0.5*std::pow(-myE2 + yE2, 2)/std::pow(sigma_yE2, 2);
            break;
        case 165:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE2, 2)) + 0.5*std::pow(-mypE2 + ypE2, 2)/std::pow(sigma_ypE2, 2);
            break;
        case 166:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yE3, 2)) + 0.5*std::pow(-myE3 + yE3, 2)/std::pow(sigma_yE3, 2);
            break;
        case 167:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yE4, 2)) + 0.5*std::pow(-myE4 + yE4, 2)/std::pow(sigma_yE4, 2);
            break;
        case 168:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE4, 2)) + 0.5*std::pow(-mypE4 + ypE4, 2)/std::pow(sigma_ypE4, 2);
            break;
        case 169:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yEv3, 2)) + 0.5*std::pow(-myEv3 + yEv3, 2)/std::pow(sigma_yEv3, 2);
            break;
        case 170:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMet, 2)) + 0.5*std::pow(-myMet + yMet, 2)/std::pow(sigma_yMet, 2);
            break;
        case 171:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yPr, 2)) + 0.5*std::pow(-myPr + yPr, 2)/std::pow(sigma_yPr, 2);
            break;
        case 172:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yFr, 2)) + 0.5*std::pow(-myFr + yFr, 2)/std::pow(sigma_yFr, 2);
            break;
        case 173:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yIr, 2)) + 0.5*std::pow(-myIr + yIr, 2)/std::pow(sigma_yIr, 2);
            break;
        case 174:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yIsr, 2)) + 0.5*std::pow(-myIsr + yIsr, 2)/std::pow(sigma_yIsr, 2);
            break;
        case 175:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yE1E1, 2)) + 0.5*std::pow(-myE1E1 + yE1E1, 2)/std::pow(sigma_yE1E1, 2);
            break;
        case 176:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yE1E2, 2)) + 0.5*std::pow(-myE1E2 + yE1E2, 2)/std::pow(sigma_yE1E2, 2);
            break;
        case 177:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yE1E3, 2)) + 0.5*std::pow(-myE1E3 + yE1E3, 2)/std::pow(sigma_yE1E3, 2);
            break;
        case 178:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yE1E4, 2)) + 0.5*std::pow(-myE1E4 + yE1E4, 2)/std::pow(sigma_yE1E4, 2);
            break;
        case 179:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yE2E2, 2)) + 0.5*std::pow(-myE2E2 + yE2E2, 2)/std::pow(sigma_yE2E2, 2);
            break;
        case 180:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yE2E3, 2)) + 0.5*std::pow(-myE2E3 + yE2E3, 2)/std::pow(sigma_yE2E3, 2);
            break;
        case 181:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yE2E4, 2)) + 0.5*std::pow(-myE2E4 + yE2E4, 2)/std::pow(sigma_yE2E4, 2);
            break;
        case 182:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yE3E4, 2)) + 0.5*std::pow(-myE3E4 + yE3E4, 2)/std::pow(sigma_yE3E4, 2);
            break;
        case 183:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yE4E4, 2)) + 0.5*std::pow(-myE4E4 + yE4E4, 2)/std::pow(sigma_yE4E4, 2);
            break;
        case 184:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMet_Met, 2)) + 0.5*std::pow(-myMet_Met + yMet_Met, 2)/std::pow(sigma_yMet_Met, 2);
            break;
        case 185:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yFrFr, 2)) + 0.5*std::pow(-myFrFr + yFrFr, 2)/std::pow(sigma_yFrFr, 2);
            break;
        case 186:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yIrIr, 2)) + 0.5*std::pow(-myIrIr + yIrIr, 2)/std::pow(sigma_yIrIr, 2);
            break;
        case 187:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yIsr_Isr, 2)) + 0.5*std::pow(-myIsr_Isr + yIsr_Isr, 2)/std::pow(sigma_yIsr_Isr, 2);
            break;
        case 188:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yEE1, 2)) + 0.5*std::pow(-myEE1 + yEE1, 2)/std::pow(sigma_yEE1, 2);
            break;
        case 189:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yHE3, 2)) + 0.5*std::pow(-myHE3 + yHE3, 2)/std::pow(sigma_yHE3, 2);
            break;
        case 190:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yHE4, 2)) + 0.5*std::pow(-myHE4 + yHE4, 2)/std::pow(sigma_yHE4, 2);
            break;
        case 191:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yHGF_Met, 2)) + 0.5*std::pow(-myHGF_Met + yHGF_Met, 2)/std::pow(sigma_yHGF_Met, 2);
            break;
        case 192:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yPPr, 2)) + 0.5*std::pow(-myPPr + yPPr, 2)/std::pow(sigma_yPPr, 2);
            break;
        case 193:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yFFr, 2)) + 0.5*std::pow(-myFFr + yFFr, 2)/std::pow(sigma_yFFr, 2);
            break;
        case 194:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yEE1E2, 2)) + 0.5*std::pow(-myEE1E2 + yEE1E2, 2)/std::pow(sigma_yEE1E2, 2);
            break;
        case 195:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yEE1Ev3, 2)) + 0.5*std::pow(-myEE1Ev3 + yEE1Ev3, 2)/std::pow(sigma_yEE1Ev3, 2);
            break;
        case 196:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yEE1E1, 2)) + 0.5*std::pow(-myEE1E1 + yEE1E1, 2)/std::pow(sigma_yEE1E1, 2);
            break;
        case 197:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yEE1E3, 2)) + 0.5*std::pow(-myEE1E3 + yEE1E3, 2)/std::pow(sigma_yEE1E3, 2);
            break;
        case 198:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yEE1E4, 2)) + 0.5*std::pow(-myEE1E4 + yEE1E4, 2)/std::pow(sigma_yEE1E4, 2);
            break;
        case 199:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yE2HE3, 2)) + 0.5*std::pow(-myE2HE3 + yE2HE3, 2)/std::pow(sigma_yE2HE3, 2);
            break;
        case 200:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yE1HE3, 2)) + 0.5*std::pow(-myE1HE3 + yE1HE3, 2)/std::pow(sigma_yE1HE3, 2);
            break;
        case 201:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yHE3E3, 2)) + 0.5*std::pow(-myHE3E3 + yHE3E3, 2)/std::pow(sigma_yHE3E3, 2);
            break;
        case 202:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yHE3Ev3, 2)) + 0.5*std::pow(-myHE3Ev3 + yHE3Ev3, 2)/std::pow(sigma_yHE3Ev3, 2);
            break;
        case 203:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yHE3E4, 2)) + 0.5*std::pow(-myHE3E4 + yHE3E4, 2)/std::pow(sigma_yHE3E4, 2);
            break;
        case 204:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yE2HE4, 2)) + 0.5*std::pow(-myE2HE4 + yE2HE4, 2)/std::pow(sigma_yE2HE4, 2);
            break;
        case 205:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yHE4Ev3, 2)) + 0.5*std::pow(-myHE4Ev3 + yHE4Ev3, 2)/std::pow(sigma_yHE4Ev3, 2);
            break;
        case 206:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yE1HE4, 2)) + 0.5*std::pow(-myE1HE4 + yE1HE4, 2)/std::pow(sigma_yE1HE4, 2);
            break;
        case 207:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yE3HE4, 2)) + 0.5*std::pow(-myE3HE4 + yE3HE4, 2)/std::pow(sigma_yE3HE4, 2);
            break;
        case 208:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yHE4E4, 2)) + 0.5*std::pow(-myHE4E4 + yHE4E4, 2)/std::pow(sigma_yHE4E4, 2);
            break;
        case 209:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yHGF_Met_Met, 2)) + 0.5*std::pow(-myHGF_Met_Met + yHGF_Met_Met, 2)/std::pow(sigma_yHGF_Met_Met, 2);
            break;
        case 210:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yPPrPr, 2)) + 0.5*std::pow(-myPPrPr + yPPrPr, 2)/std::pow(sigma_yPPrPr, 2);
            break;
        case 211:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yFFrFr, 2)) + 0.5*std::pow(-myFFrFr + yFFrFr, 2)/std::pow(sigma_yFFrFr, 2);
            break;
        case 212:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yIIrIr, 2)) + 0.5*std::pow(-myIIrIr + yIIrIr, 2)/std::pow(sigma_yIIrIr, 2);
            break;
        case 213:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yINS_Isr_Isr, 2)) + 0.5*std::pow(-myINS_Isr_Isr + yINS_Isr_Isr, 2)/std::pow(sigma_yINS_Isr_Isr, 2);
            break;
        case 214:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yEE1EE1, 2)) + 0.5*std::pow(-myEE1EE1 + yEE1EE1, 2)/std::pow(sigma_yEE1EE1, 2);
            break;
        case 215:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yEE1HE3, 2)) + 0.5*std::pow(-myEE1HE3 + yEE1HE3, 2)/std::pow(sigma_yEE1HE3, 2);
            break;
        case 216:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yEE1HE4, 2)) + 0.5*std::pow(-myEE1HE4 + yEE1HE4, 2)/std::pow(sigma_yEE1HE4, 2);
            break;
        case 217:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yHE3HE3, 2)) + 0.5*std::pow(-myHE3HE3 + yHE3HE3, 2)/std::pow(sigma_yHE3HE3, 2);
            break;
        case 218:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yHE3HE4, 2)) + 0.5*std::pow(-myHE3HE4 + yHE3HE4, 2)/std::pow(sigma_yHE3HE4, 2);
            break;
        case 219:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yHE4HE4, 2)) + 0.5*std::pow(-myHE4HE4 + yHE4HE4, 2)/std::pow(sigma_yHE4HE4, 2);
            break;
        case 220:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yHGF_Met_HGF_Met, 2)) + 0.5*std::pow(-myHGF_Met_HGF_Met + yHGF_Met_HGF_Met, 2)/std::pow(sigma_yHGF_Met_HGF_Met, 2);
            break;
        case 221:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yPPrPPr, 2)) + 0.5*std::pow(-myPPrPPr + yPPrPPr, 2)/std::pow(sigma_yPPrPPr, 2);
            break;
        case 222:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yFFrFFr, 2)) + 0.5*std::pow(-myFFrFFr + yFFrFFr, 2)/std::pow(sigma_yFFrFFr, 2);
            break;
        case 223:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yIIrIrI, 2)) + 0.5*std::pow(-myIIrIrI + yIIrIrI, 2)/std::pow(sigma_yIIrIrI, 2);
            break;
        case 224:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yINS_Isr_Isr_INS, 2)) + 0.5*std::pow(-myINS_Isr_Isr_INS + yINS_Isr_Isr_INS, 2)/std::pow(sigma_yINS_Isr_Isr_INS, 2);
            break;
        case 225:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yE1_ppERK, 2)) + 0.5*std::pow(-myE1_ppERK + yE1_ppERK, 2)/std::pow(sigma_yE1_ppERK, 2);
            break;
        case 226:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yE2_ppERK, 2)) + 0.5*std::pow(-myE2_ppERK + yE2_ppERK, 2)/std::pow(sigma_yE2_ppERK, 2);
            break;
        case 227:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yE4_ppERK, 2)) + 0.5*std::pow(-myE4_ppERK + yE4_ppERK, 2)/std::pow(sigma_yE4_ppERK, 2);
            break;
        case 228:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E2, 2)) + 0.5*std::pow(-mypEE1E2 + ypEE1E2, 2)/std::pow(sigma_ypEE1E2, 2);
            break;
        case 229:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1Ev3, 2)) + 0.5*std::pow(-mypEE1Ev3 + ypEE1Ev3, 2)/std::pow(sigma_ypEE1Ev3, 2);
            break;
        case 230:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E1, 2)) + 0.5*std::pow(-mypEE1E1 + ypEE1E1, 2)/std::pow(sigma_ypEE1E1, 2);
            break;
        case 231:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1EE1, 2)) + 0.5*std::pow(-mypEE1EE1 + ypEE1EE1, 2)/std::pow(sigma_ypEE1EE1, 2);
            break;
        case 232:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E3, 2)) + 0.5*std::pow(-mypEE1E3 + ypEE1E3, 2)/std::pow(sigma_ypEE1E3, 2);
            break;
        case 233:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1HE3, 2)) + 0.5*std::pow(-mypEE1HE3 + ypEE1HE3, 2)/std::pow(sigma_ypEE1HE3, 2);
            break;
        case 234:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E4, 2)) + 0.5*std::pow(-mypEE1E4 + ypEE1E4, 2)/std::pow(sigma_ypEE1E4, 2);
            break;
        case 235:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1HE4, 2)) + 0.5*std::pow(-mypEE1HE4 + ypEE1HE4, 2)/std::pow(sigma_ypEE1HE4, 2);
            break;
        case 236:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE2HE3, 2)) + 0.5*std::pow(-mypE2HE3 + ypE2HE3, 2)/std::pow(sigma_ypE2HE3, 2);
            break;
        case 237:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE3Ev3, 2)) + 0.5*std::pow(-mypHE3Ev3 + ypHE3Ev3, 2)/std::pow(sigma_ypHE3Ev3, 2);
            break;
        case 238:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE1HE3, 2)) + 0.5*std::pow(-mypE1HE3 + ypE1HE3, 2)/std::pow(sigma_ypE1HE3, 2);
            break;
        case 239:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE3E4, 2)) + 0.5*std::pow(-mypHE3E4 + ypHE3E4, 2)/std::pow(sigma_ypHE3E4, 2);
            break;
        case 240:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE3HE4, 2)) + 0.5*std::pow(-mypHE3HE4 + ypHE3HE4, 2)/std::pow(sigma_ypHE3HE4, 2);
            break;
        case 241:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE2HE4, 2)) + 0.5*std::pow(-mypE2HE4 + ypE2HE4, 2)/std::pow(sigma_ypE2HE4, 2);
            break;
        case 242:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE4Ev3, 2)) + 0.5*std::pow(-mypHE4Ev3 + ypHE4Ev3, 2)/std::pow(sigma_ypHE4Ev3, 2);
            break;
        case 243:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE1HE4, 2)) + 0.5*std::pow(-mypE1HE4 + ypE1HE4, 2)/std::pow(sigma_ypE1HE4, 2);
            break;
        case 244:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE3HE4, 2)) + 0.5*std::pow(-mypE3HE4 + ypE3HE4, 2)/std::pow(sigma_ypE3HE4, 2);
            break;
        case 245:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE4E4, 2)) + 0.5*std::pow(-mypHE4E4 + ypHE4E4, 2)/std::pow(sigma_ypHE4E4, 2);
            break;
        case 246:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE4HE4, 2)) + 0.5*std::pow(-mypHE4HE4 + ypHE4HE4, 2)/std::pow(sigma_ypHE4HE4, 2);
            break;
        case 247:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHGF_Met_Met, 2)) + 0.5*std::pow(-mypHGF_Met_Met + ypHGF_Met_Met, 2)/std::pow(sigma_ypHGF_Met_Met, 2);
            break;
        case 248:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHGF_Met_HGF_Met, 2)) + 0.5*std::pow(-mypHGF_Met_HGF_Met + ypHGF_Met_HGF_Met, 2)/std::pow(sigma_ypHGF_Met_HGF_Met, 2);
            break;
        case 249:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypPPrPPr, 2)) + 0.5*std::pow(-mypPPrPPr + ypPPrPPr, 2)/std::pow(sigma_ypPPrPPr, 2);
            break;
        case 250:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypPPrPr, 2)) + 0.5*std::pow(-mypPPrPr + ypPPrPr, 2)/std::pow(sigma_ypPPrPr, 2);
            break;
        case 251:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypFFrFFr, 2)) + 0.5*std::pow(-mypFFrFFr + ypFFrFFr, 2)/std::pow(sigma_ypFFrFFr, 2);
            break;
        case 252:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypFFrFr, 2)) + 0.5*std::pow(-mypFFrFr + ypFFrFr, 2)/std::pow(sigma_ypFFrFr, 2);
            break;
        case 253:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypIIrIr, 2)) + 0.5*std::pow(-mypIIrIr + ypIIrIr, 2)/std::pow(sigma_ypIIrIr, 2);
            break;
        case 254:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypINS_Isr_Isr, 2)) + 0.5*std::pow(-mypINS_Isr_Isr + ypINS_Isr_Isr, 2)/std::pow(sigma_ypINS_Isr_Isr, 2);
            break;
        case 255:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypIIrIrI, 2)) + 0.5*std::pow(-mypIIrIrI + ypIIrIrI, 2)/std::pow(sigma_ypIIrIrI, 2);
            break;
        case 256:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypINS_Isr_Isr_INS, 2)) + 0.5*std::pow(-mypINS_Isr_Isr_INS + ypINS_Isr_Isr_INS, 2)/std::pow(sigma_ypINS_Isr_Isr_INS, 2);
            break;
        case 257:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypIIrIr_IRS, 2)) + 0.5*std::pow(-mypIIrIr_IRS + ypIIrIr_IRS, 2)/std::pow(sigma_ypIIrIr_IRS, 2);
            break;
        case 258:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypINS_Isr_Isr_IRS, 2)) + 0.5*std::pow(-mypINS_Isr_Isr_IRS + ypINS_Isr_Isr_IRS, 2)/std::pow(sigma_ypINS_Isr_Isr_IRS, 2);
            break;
        case 259:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypIIrIrI_IRS, 2)) + 0.5*std::pow(-mypIIrIrI_IRS + ypIIrIrI_IRS, 2)/std::pow(sigma_ypIIrIrI_IRS, 2);
            break;
        case 260:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypINS_Isr_Isr_INS_IRS, 2)) + 0.5*std::pow(-mypINS_Isr_Isr_INS_IRS + ypINS_Isr_Isr_INS_IRS, 2)/std::pow(sigma_ypINS_Isr_Isr_INS_IRS, 2);
            break;
        case 261:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ySp_EE1E2, 2)) + 0.5*std::pow(-mySp_EE1E2 + ySp_EE1E2, 2)/std::pow(sigma_ySp_EE1E2, 2);
            break;
        case 262:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ySp_EE1Ev3, 2)) + 0.5*std::pow(-mySp_EE1Ev3 + ySp_EE1Ev3, 2)/std::pow(sigma_ySp_EE1Ev3, 2);
            break;
        case 263:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ySp_EE1E1, 2)) + 0.5*std::pow(-mySp_EE1E1 + ySp_EE1E1, 2)/std::pow(sigma_ySp_EE1E1, 2);
            break;
        case 264:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ySp_EE1EE1, 2)) + 0.5*std::pow(-mySp_EE1EE1 + ySp_EE1EE1, 2)/std::pow(sigma_ySp_EE1EE1, 2);
            break;
        case 265:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ySp_EE1E3, 2)) + 0.5*std::pow(-mySp_EE1E3 + ySp_EE1E3, 2)/std::pow(sigma_ySp_EE1E3, 2);
            break;
        case 266:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ySp_EE1HE3, 2)) + 0.5*std::pow(-mySp_EE1HE3 + ySp_EE1HE3, 2)/std::pow(sigma_ySp_EE1HE3, 2);
            break;
        case 267:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ySp_EE1E4, 2)) + 0.5*std::pow(-mySp_EE1E4 + ySp_EE1E4, 2)/std::pow(sigma_ySp_EE1E4, 2);
            break;
        case 268:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ySp_EE1HE4, 2)) + 0.5*std::pow(-mySp_EE1HE4 + ySp_EE1HE4, 2)/std::pow(sigma_ySp_EE1HE4, 2);
            break;
        case 269:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ySp_E2HE3, 2)) + 0.5*std::pow(-mySp_E2HE3 + ySp_E2HE3, 2)/std::pow(sigma_ySp_E2HE3, 2);
            break;
        case 270:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ySp_HE3Ev3, 2)) + 0.5*std::pow(-mySp_HE3Ev3 + ySp_HE3Ev3, 2)/std::pow(sigma_ySp_HE3Ev3, 2);
            break;
        case 271:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ySp_E1HE3, 2)) + 0.5*std::pow(-mySp_E1HE3 + ySp_E1HE3, 2)/std::pow(sigma_ySp_E1HE3, 2);
            break;
        case 272:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ySp_HE3E4, 2)) + 0.5*std::pow(-mySp_HE3E4 + ySp_HE3E4, 2)/std::pow(sigma_ySp_HE3E4, 2);
            break;
        case 273:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ySp_HE3HE4, 2)) + 0.5*std::pow(-mySp_HE3HE4 + ySp_HE3HE4, 2)/std::pow(sigma_ySp_HE3HE4, 2);
            break;
        case 274:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ySp_E2HE4, 2)) + 0.5*std::pow(-mySp_E2HE4 + ySp_E2HE4, 2)/std::pow(sigma_ySp_E2HE4, 2);
            break;
        case 275:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ySp_HE4Ev3, 2)) + 0.5*std::pow(-mySp_HE4Ev3 + ySp_HE4Ev3, 2)/std::pow(sigma_ySp_HE4Ev3, 2);
            break;
        case 276:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ySp_E1HE4, 2)) + 0.5*std::pow(-mySp_E1HE4 + ySp_E1HE4, 2)/std::pow(sigma_ySp_E1HE4, 2);
            break;
        case 277:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ySp_E3HE4, 2)) + 0.5*std::pow(-mySp_E3HE4 + ySp_E3HE4, 2)/std::pow(sigma_ySp_E3HE4, 2);
            break;
        case 278:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ySp_HE4E4, 2)) + 0.5*std::pow(-mySp_HE4E4 + ySp_HE4E4, 2)/std::pow(sigma_ySp_HE4E4, 2);
            break;
        case 279:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ySp_HE4HE4, 2)) + 0.5*std::pow(-mySp_HE4HE4 + ySp_HE4HE4, 2)/std::pow(sigma_ySp_HE4HE4, 2);
            break;
        case 280:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ySp_HGF_Met_Met, 2)) + 0.5*std::pow(-mySp_HGF_Met_Met + ySp_HGF_Met_Met, 2)/std::pow(sigma_ySp_HGF_Met_Met, 2);
            break;
        case 281:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ySp_HGF_Met_HGF_Met, 2)) + 0.5*std::pow(-mySp_HGF_Met_HGF_Met + ySp_HGF_Met_HGF_Met, 2)/std::pow(sigma_ySp_HGF_Met_HGF_Met, 2);
            break;
        case 282:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ySp_PPrPPr, 2)) + 0.5*std::pow(-mySp_PPrPPr + ySp_PPrPPr, 2)/std::pow(sigma_ySp_PPrPPr, 2);
            break;
        case 283:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ySp_PPrPr, 2)) + 0.5*std::pow(-mySp_PPrPr + ySp_PPrPr, 2)/std::pow(sigma_ySp_PPrPr, 2);
            break;
        case 284:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ySp_FFrFFr, 2)) + 0.5*std::pow(-mySp_FFrFFr + ySp_FFrFFr, 2)/std::pow(sigma_ySp_FFrFFr, 2);
            break;
        case 285:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ySp_FFrFr, 2)) + 0.5*std::pow(-mySp_FFrFr + ySp_FFrFr, 2)/std::pow(sigma_ySp_FFrFr, 2);
            break;
        case 286:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ySp_IIrIr, 2)) + 0.5*std::pow(-mySp_IIrIr + ySp_IIrIr, 2)/std::pow(sigma_ySp_IIrIr, 2);
            break;
        case 287:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ySp_INS_Isr_Isr, 2)) + 0.5*std::pow(-mySp_INS_Isr_Isr + ySp_INS_Isr_Isr, 2)/std::pow(sigma_ySp_INS_Isr_Isr, 2);
            break;
        case 288:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ySp_IIrIrI, 2)) + 0.5*std::pow(-mySp_IIrIrI + ySp_IIrIrI, 2)/std::pow(sigma_ySp_IIrIrI, 2);
            break;
        case 289:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ySp_INS_Isr_Isr_INS, 2)) + 0.5*std::pow(-mySp_INS_Isr_Isr_INS + ySp_INS_Isr_Isr_INS, 2)/std::pow(sigma_ySp_INS_Isr_Isr_INS, 2);
            break;
        case 290:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yEE1E2int, 2)) + 0.5*std::pow(-myEE1E2int + yEE1E2int, 2)/std::pow(sigma_yEE1E2int, 2);
            break;
        case 291:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yEE1Ev3int, 2)) + 0.5*std::pow(-myEE1Ev3int + yEE1Ev3int, 2)/std::pow(sigma_yEE1Ev3int, 2);
            break;
        case 292:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yEE1E1int, 2)) + 0.5*std::pow(-myEE1E1int + yEE1E1int, 2)/std::pow(sigma_yEE1E1int, 2);
            break;
        case 293:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yEE1EE1int, 2)) + 0.5*std::pow(-myEE1EE1int + yEE1EE1int, 2)/std::pow(sigma_yEE1EE1int, 2);
            break;
        case 294:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yEE1E3int, 2)) + 0.5*std::pow(-myEE1E3int + yEE1E3int, 2)/std::pow(sigma_yEE1E3int, 2);
            break;
        case 295:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yEE1HE3int, 2)) + 0.5*std::pow(-myEE1HE3int + yEE1HE3int, 2)/std::pow(sigma_yEE1HE3int, 2);
            break;
        case 296:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yEE1E4int, 2)) + 0.5*std::pow(-myEE1E4int + yEE1E4int, 2)/std::pow(sigma_yEE1E4int, 2);
            break;
        case 297:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yEE1HE4int, 2)) + 0.5*std::pow(-myEE1HE4int + yEE1HE4int, 2)/std::pow(sigma_yEE1HE4int, 2);
            break;
        case 298:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yE2HE3int, 2)) + 0.5*std::pow(-myE2HE3int + yE2HE3int, 2)/std::pow(sigma_yE2HE3int, 2);
            break;
        case 299:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yHE3Ev3int, 2)) + 0.5*std::pow(-myHE3Ev3int + yHE3Ev3int, 2)/std::pow(sigma_yHE3Ev3int, 2);
            break;
        case 300:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yE1HE3int, 2)) + 0.5*std::pow(-myE1HE3int + yE1HE3int, 2)/std::pow(sigma_yE1HE3int, 2);
            break;
        case 301:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yHE3E4int, 2)) + 0.5*std::pow(-myHE3E4int + yHE3E4int, 2)/std::pow(sigma_yHE3E4int, 2);
            break;
        case 302:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yHE3HE4int, 2)) + 0.5*std::pow(-myHE3HE4int + yHE3HE4int, 2)/std::pow(sigma_yHE3HE4int, 2);
            break;
        case 303:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yE2HE4int, 2)) + 0.5*std::pow(-myE2HE4int + yE2HE4int, 2)/std::pow(sigma_yE2HE4int, 2);
            break;
        case 304:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yHE4Ev3int, 2)) + 0.5*std::pow(-myHE4Ev3int + yHE4Ev3int, 2)/std::pow(sigma_yHE4Ev3int, 2);
            break;
        case 305:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yE1HE4int, 2)) + 0.5*std::pow(-myE1HE4int + yE1HE4int, 2)/std::pow(sigma_yE1HE4int, 2);
            break;
        case 306:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yE3HE4int, 2)) + 0.5*std::pow(-myE3HE4int + yE3HE4int, 2)/std::pow(sigma_yE3HE4int, 2);
            break;
        case 307:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yHE4E4int, 2)) + 0.5*std::pow(-myHE4E4int + yHE4E4int, 2)/std::pow(sigma_yHE4E4int, 2);
            break;
        case 308:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yHE4HE4int, 2)) + 0.5*std::pow(-myHE4HE4int + yHE4HE4int, 2)/std::pow(sigma_yHE4HE4int, 2);
            break;
        case 309:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yHGF_Met_Metint, 2)) + 0.5*std::pow(-myHGF_Met_Metint + yHGF_Met_Metint, 2)/std::pow(sigma_yHGF_Met_Metint, 2);
            break;
        case 310:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yHGF_Met_HGF_Metint, 2)) + 0.5*std::pow(-myHGF_Met_HGF_Metint + yHGF_Met_HGF_Metint, 2)/std::pow(sigma_yHGF_Met_HGF_Metint, 2);
            break;
        case 311:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yPPrPPrint, 2)) + 0.5*std::pow(-myPPrPPrint + yPPrPPrint, 2)/std::pow(sigma_yPPrPPrint, 2);
            break;
        case 312:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yPPrPrint, 2)) + 0.5*std::pow(-myPPrPrint + yPPrPrint, 2)/std::pow(sigma_yPPrPrint, 2);
            break;
        case 313:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yFFrFFrint, 2)) + 0.5*std::pow(-myFFrFFrint + yFFrFFrint, 2)/std::pow(sigma_yFFrFFrint, 2);
            break;
        case 314:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yFFrFrint, 2)) + 0.5*std::pow(-myFFrFrint + yFFrFrint, 2)/std::pow(sigma_yFFrFrint, 2);
            break;
        case 315:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yIIrIr_int, 2)) + 0.5*std::pow(-myIIrIr_int + yIIrIr_int, 2)/std::pow(sigma_yIIrIr_int, 2);
            break;
        case 316:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yINS_Isr_Isr_int, 2)) + 0.5*std::pow(-myINS_Isr_Isr_int + yINS_Isr_Isr_int, 2)/std::pow(sigma_yINS_Isr_Isr_int, 2);
            break;
        case 317:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yIIrIrI_int, 2)) + 0.5*std::pow(-myIIrIrI_int + yIIrIrI_int, 2)/std::pow(sigma_yIIrIrI_int, 2);
            break;
        case 318:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yINS_Isr_Isr_INS_int, 2)) + 0.5*std::pow(-myINS_Isr_Isr_INS_int + yINS_Isr_Isr_INS_int, 2)/std::pow(sigma_yINS_Isr_Isr_INS_int, 2);
            break;
        case 319:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E2int, 2)) + 0.5*std::pow(-mypEE1E2int + ypEE1E2int, 2)/std::pow(sigma_ypEE1E2int, 2);
            break;
        case 320:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1Ev3int, 2)) + 0.5*std::pow(-mypEE1Ev3int + ypEE1Ev3int, 2)/std::pow(sigma_ypEE1Ev3int, 2);
            break;
        case 321:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E1int, 2)) + 0.5*std::pow(-mypEE1E1int + ypEE1E1int, 2)/std::pow(sigma_ypEE1E1int, 2);
            break;
        case 322:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1EE1int, 2)) + 0.5*std::pow(-mypEE1EE1int + ypEE1EE1int, 2)/std::pow(sigma_ypEE1EE1int, 2);
            break;
        case 323:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E3int, 2)) + 0.5*std::pow(-mypEE1E3int + ypEE1E3int, 2)/std::pow(sigma_ypEE1E3int, 2);
            break;
        case 324:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1HE3int, 2)) + 0.5*std::pow(-mypEE1HE3int + ypEE1HE3int, 2)/std::pow(sigma_ypEE1HE3int, 2);
            break;
        case 325:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E4int, 2)) + 0.5*std::pow(-mypEE1E4int + ypEE1E4int, 2)/std::pow(sigma_ypEE1E4int, 2);
            break;
        case 326:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1HE4int, 2)) + 0.5*std::pow(-mypEE1HE4int + ypEE1HE4int, 2)/std::pow(sigma_ypEE1HE4int, 2);
            break;
        case 327:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE2HE3int, 2)) + 0.5*std::pow(-mypE2HE3int + ypE2HE3int, 2)/std::pow(sigma_ypE2HE3int, 2);
            break;
        case 328:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE3Ev3int, 2)) + 0.5*std::pow(-mypHE3Ev3int + ypHE3Ev3int, 2)/std::pow(sigma_ypHE3Ev3int, 2);
            break;
        case 329:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE1HE3int, 2)) + 0.5*std::pow(-mypE1HE3int + ypE1HE3int, 2)/std::pow(sigma_ypE1HE3int, 2);
            break;
        case 330:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE3E4int, 2)) + 0.5*std::pow(-mypHE3E4int + ypHE3E4int, 2)/std::pow(sigma_ypHE3E4int, 2);
            break;
        case 331:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE3HE4int, 2)) + 0.5*std::pow(-mypHE3HE4int + ypHE3HE4int, 2)/std::pow(sigma_ypHE3HE4int, 2);
            break;
        case 332:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE2HE4int, 2)) + 0.5*std::pow(-mypE2HE4int + ypE2HE4int, 2)/std::pow(sigma_ypE2HE4int, 2);
            break;
        case 333:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE4Ev3int, 2)) + 0.5*std::pow(-mypHE4Ev3int + ypHE4Ev3int, 2)/std::pow(sigma_ypHE4Ev3int, 2);
            break;
        case 334:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE1HE4int, 2)) + 0.5*std::pow(-mypE1HE4int + ypE1HE4int, 2)/std::pow(sigma_ypE1HE4int, 2);
            break;
        case 335:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE3HE4int, 2)) + 0.5*std::pow(-mypE3HE4int + ypE3HE4int, 2)/std::pow(sigma_ypE3HE4int, 2);
            break;
        case 336:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE4E4int, 2)) + 0.5*std::pow(-mypHE4E4int + ypHE4E4int, 2)/std::pow(sigma_ypHE4E4int, 2);
            break;
        case 337:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE4HE4int, 2)) + 0.5*std::pow(-mypHE4HE4int + ypHE4HE4int, 2)/std::pow(sigma_ypHE4HE4int, 2);
            break;
        case 338:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHGF_Met_Metint, 2)) + 0.5*std::pow(-mypHGF_Met_Metint + ypHGF_Met_Metint, 2)/std::pow(sigma_ypHGF_Met_Metint, 2);
            break;
        case 339:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHGF_Met_HGF_Metint, 2)) + 0.5*std::pow(-mypHGF_Met_HGF_Metint + ypHGF_Met_HGF_Metint, 2)/std::pow(sigma_ypHGF_Met_HGF_Metint, 2);
            break;
        case 340:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypPPrPPrint, 2)) + 0.5*std::pow(-mypPPrPPrint + ypPPrPPrint, 2)/std::pow(sigma_ypPPrPPrint, 2);
            break;
        case 341:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypPPrPrint, 2)) + 0.5*std::pow(-mypPPrPrint + ypPPrPrint, 2)/std::pow(sigma_ypPPrPrint, 2);
            break;
        case 342:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypFFrFFrint, 2)) + 0.5*std::pow(-mypFFrFFrint + ypFFrFFrint, 2)/std::pow(sigma_ypFFrFFrint, 2);
            break;
        case 343:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypFFrFrint, 2)) + 0.5*std::pow(-mypFFrFrint + ypFFrFrint, 2)/std::pow(sigma_ypFFrFrint, 2);
            break;
        case 344:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypIIrIr_int, 2)) + 0.5*std::pow(-mypIIrIr_int + ypIIrIr_int, 2)/std::pow(sigma_ypIIrIr_int, 2);
            break;
        case 345:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypINS_Isr_Isr_int, 2)) + 0.5*std::pow(-mypINS_Isr_Isr_int + ypINS_Isr_Isr_int, 2)/std::pow(sigma_ypINS_Isr_Isr_int, 2);
            break;
        case 346:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypIIrIrI_int, 2)) + 0.5*std::pow(-mypIIrIrI_int + ypIIrIrI_int, 2)/std::pow(sigma_ypIIrIrI_int, 2);
            break;
        case 347:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypINS_Isr_Isr_INS_int, 2)) + 0.5*std::pow(-mypINS_Isr_Isr_INS_int + ypINS_Isr_Isr_INS_int, 2)/std::pow(sigma_ypINS_Isr_Isr_INS_int, 2);
            break;
        case 348:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypIIrIr_int_IRS, 2)) + 0.5*std::pow(-mypIIrIr_int_IRS + ypIIrIr_int_IRS, 2)/std::pow(sigma_ypIIrIr_int_IRS, 2);
            break;
        case 349:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypINS_Isr_Isr_int_IRS, 2)) + 0.5*std::pow(-mypINS_Isr_Isr_int_IRS + ypINS_Isr_Isr_int_IRS, 2)/std::pow(sigma_ypINS_Isr_Isr_int_IRS, 2);
            break;
        case 350:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypIIrIrI_int_IRS, 2)) + 0.5*std::pow(-mypIIrIrI_int_IRS + ypIIrIrI_int_IRS, 2)/std::pow(sigma_ypIIrIrI_int_IRS, 2);
            break;
        case 351:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypINS_Isr_Isr_INS_int_IRS, 2)) + 0.5*std::pow(-mypINS_Isr_Isr_INS_int_IRS + ypINS_Isr_Isr_INS_int_IRS, 2)/std::pow(sigma_ypINS_Isr_Isr_INS_int_IRS, 2);
            break;
        case 352:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E2_G2_SOS, 2)) + 0.5*std::pow(-mypEE1E2_G2_SOS + ypEE1E2_G2_SOS, 2)/std::pow(sigma_ypEE1E2_G2_SOS, 2);
            break;
        case 353:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1Ev3_G2_SOS, 2)) + 0.5*std::pow(-mypEE1Ev3_G2_SOS + ypEE1Ev3_G2_SOS, 2)/std::pow(sigma_ypEE1Ev3_G2_SOS, 2);
            break;
        case 354:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E1_G2_SOS, 2)) + 0.5*std::pow(-mypEE1E1_G2_SOS + ypEE1E1_G2_SOS, 2)/std::pow(sigma_ypEE1E1_G2_SOS, 2);
            break;
        case 355:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1EE1_G2_SOS, 2)) + 0.5*std::pow(-mypEE1EE1_G2_SOS + ypEE1EE1_G2_SOS, 2)/std::pow(sigma_ypEE1EE1_G2_SOS, 2);
            break;
        case 356:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E3_G2_SOS, 2)) + 0.5*std::pow(-mypEE1E3_G2_SOS + ypEE1E3_G2_SOS, 2)/std::pow(sigma_ypEE1E3_G2_SOS, 2);
            break;
        case 357:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1HE3_G2_SOS, 2)) + 0.5*std::pow(-mypEE1HE3_G2_SOS + ypEE1HE3_G2_SOS, 2)/std::pow(sigma_ypEE1HE3_G2_SOS, 2);
            break;
        case 358:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E4_G2_SOS, 2)) + 0.5*std::pow(-mypEE1E4_G2_SOS + ypEE1E4_G2_SOS, 2)/std::pow(sigma_ypEE1E4_G2_SOS, 2);
            break;
        case 359:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1HE4_G2_SOS, 2)) + 0.5*std::pow(-mypEE1HE4_G2_SOS + ypEE1HE4_G2_SOS, 2)/std::pow(sigma_ypEE1HE4_G2_SOS, 2);
            break;
        case 360:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE2HE3_G2_SOS, 2)) + 0.5*std::pow(-mypE2HE3_G2_SOS + ypE2HE3_G2_SOS, 2)/std::pow(sigma_ypE2HE3_G2_SOS, 2);
            break;
        case 361:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE3Ev3_G2_SOS, 2)) + 0.5*std::pow(-mypHE3Ev3_G2_SOS + ypHE3Ev3_G2_SOS, 2)/std::pow(sigma_ypHE3Ev3_G2_SOS, 2);
            break;
        case 362:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE1HE3_G2_SOS, 2)) + 0.5*std::pow(-mypE1HE3_G2_SOS + ypE1HE3_G2_SOS, 2)/std::pow(sigma_ypE1HE3_G2_SOS, 2);
            break;
        case 363:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE3E4_G2_SOS, 2)) + 0.5*std::pow(-mypHE3E4_G2_SOS + ypHE3E4_G2_SOS, 2)/std::pow(sigma_ypHE3E4_G2_SOS, 2);
            break;
        case 364:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE3HE4_G2_SOS, 2)) + 0.5*std::pow(-mypHE3HE4_G2_SOS + ypHE3HE4_G2_SOS, 2)/std::pow(sigma_ypHE3HE4_G2_SOS, 2);
            break;
        case 365:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE2HE4_G2_SOS, 2)) + 0.5*std::pow(-mypE2HE4_G2_SOS + ypE2HE4_G2_SOS, 2)/std::pow(sigma_ypE2HE4_G2_SOS, 2);
            break;
        case 366:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE4Ev3_G2_SOS, 2)) + 0.5*std::pow(-mypHE4Ev3_G2_SOS + ypHE4Ev3_G2_SOS, 2)/std::pow(sigma_ypHE4Ev3_G2_SOS, 2);
            break;
        case 367:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE1HE4_G2_SOS, 2)) + 0.5*std::pow(-mypE1HE4_G2_SOS + ypE1HE4_G2_SOS, 2)/std::pow(sigma_ypE1HE4_G2_SOS, 2);
            break;
        case 368:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE3HE4_G2_SOS, 2)) + 0.5*std::pow(-mypE3HE4_G2_SOS + ypE3HE4_G2_SOS, 2)/std::pow(sigma_ypE3HE4_G2_SOS, 2);
            break;
        case 369:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE4E4_G2_SOS, 2)) + 0.5*std::pow(-mypHE4E4_G2_SOS + ypHE4E4_G2_SOS, 2)/std::pow(sigma_ypHE4E4_G2_SOS, 2);
            break;
        case 370:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE4HE4_G2_SOS, 2)) + 0.5*std::pow(-mypHE4HE4_G2_SOS + ypHE4HE4_G2_SOS, 2)/std::pow(sigma_ypHE4HE4_G2_SOS, 2);
            break;
        case 371:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHGF_Met_Met_G2_SOS, 2)) + 0.5*std::pow(-mypHGF_Met_Met_G2_SOS + ypHGF_Met_Met_G2_SOS, 2)/std::pow(sigma_ypHGF_Met_Met_G2_SOS, 2);
            break;
        case 372:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHGF_Met_HGF_Met_G2_SOS, 2)) + 0.5*std::pow(-mypHGF_Met_HGF_Met_G2_SOS + ypHGF_Met_HGF_Met_G2_SOS, 2)/std::pow(sigma_ypHGF_Met_HGF_Met_G2_SOS, 2);
            break;
        case 373:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypPPrPPr_G2_SOS, 2)) + 0.5*std::pow(-mypPPrPPr_G2_SOS + ypPPrPPr_G2_SOS, 2)/std::pow(sigma_ypPPrPPr_G2_SOS, 2);
            break;
        case 374:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypPPrPr_G2_SOS, 2)) + 0.5*std::pow(-mypPPrPr_G2_SOS + ypPPrPr_G2_SOS, 2)/std::pow(sigma_ypPPrPr_G2_SOS, 2);
            break;
        case 375:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypFFrFFr_G2_SOS, 2)) + 0.5*std::pow(-mypFFrFFr_G2_SOS + ypFFrFFr_G2_SOS, 2)/std::pow(sigma_ypFFrFFr_G2_SOS, 2);
            break;
        case 376:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypFFrFr_G2_SOS, 2)) + 0.5*std::pow(-mypFFrFr_G2_SOS + ypFFrFr_G2_SOS, 2)/std::pow(sigma_ypFFrFr_G2_SOS, 2);
            break;
        case 377:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypIIrIr_IRS_G2_SOS, 2)) + 0.5*std::pow(-mypIIrIr_IRS_G2_SOS + ypIIrIr_IRS_G2_SOS, 2)/std::pow(sigma_ypIIrIr_IRS_G2_SOS, 2);
            break;
        case 378:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypINS_Isr_Isr_IRS_G2_SOS, 2)) + 0.5*std::pow(-mypINS_Isr_Isr_IRS_G2_SOS + ypINS_Isr_Isr_IRS_G2_SOS, 2)/std::pow(sigma_ypINS_Isr_Isr_IRS_G2_SOS, 2);
            break;
        case 379:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypIIrIrI_IRS_G2_SOS, 2)) + 0.5*std::pow(-mypIIrIrI_IRS_G2_SOS + ypIIrIrI_IRS_G2_SOS, 2)/std::pow(sigma_ypIIrIrI_IRS_G2_SOS, 2);
            break;
        case 380:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypINS_Isr_Isr_INS_IRS_G2_SOS, 2)) + 0.5*std::pow(-mypINS_Isr_Isr_INS_IRS_G2_SOS + ypINS_Isr_Isr_INS_IRS_G2_SOS, 2)/std::pow(sigma_ypINS_Isr_Isr_INS_IRS_G2_SOS, 2);
            break;
        case 381:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E2int_G2_SOS, 2)) + 0.5*std::pow(-mypEE1E2int_G2_SOS + ypEE1E2int_G2_SOS, 2)/std::pow(sigma_ypEE1E2int_G2_SOS, 2);
            break;
        case 382:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1Ev3int_G2_SOS, 2)) + 0.5*std::pow(-mypEE1Ev3int_G2_SOS + ypEE1Ev3int_G2_SOS, 2)/std::pow(sigma_ypEE1Ev3int_G2_SOS, 2);
            break;
        case 383:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E1int_G2_SOS, 2)) + 0.5*std::pow(-mypEE1E1int_G2_SOS + ypEE1E1int_G2_SOS, 2)/std::pow(sigma_ypEE1E1int_G2_SOS, 2);
            break;
        case 384:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1EE1int_G2_SOS, 2)) + 0.5*std::pow(-mypEE1EE1int_G2_SOS + ypEE1EE1int_G2_SOS, 2)/std::pow(sigma_ypEE1EE1int_G2_SOS, 2);
            break;
        case 385:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E3int_G2_SOS, 2)) + 0.5*std::pow(-mypEE1E3int_G2_SOS + ypEE1E3int_G2_SOS, 2)/std::pow(sigma_ypEE1E3int_G2_SOS, 2);
            break;
        case 386:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1HE3int_G2_SOS, 2)) + 0.5*std::pow(-mypEE1HE3int_G2_SOS + ypEE1HE3int_G2_SOS, 2)/std::pow(sigma_ypEE1HE3int_G2_SOS, 2);
            break;
        case 387:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E4int_G2_SOS, 2)) + 0.5*std::pow(-mypEE1E4int_G2_SOS + ypEE1E4int_G2_SOS, 2)/std::pow(sigma_ypEE1E4int_G2_SOS, 2);
            break;
        case 388:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1HE4int_G2_SOS, 2)) + 0.5*std::pow(-mypEE1HE4int_G2_SOS + ypEE1HE4int_G2_SOS, 2)/std::pow(sigma_ypEE1HE4int_G2_SOS, 2);
            break;
        case 389:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE2HE3int_G2_SOS, 2)) + 0.5*std::pow(-mypE2HE3int_G2_SOS + ypE2HE3int_G2_SOS, 2)/std::pow(sigma_ypE2HE3int_G2_SOS, 2);
            break;
        case 390:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE3Ev3int_G2_SOS, 2)) + 0.5*std::pow(-mypHE3Ev3int_G2_SOS + ypHE3Ev3int_G2_SOS, 2)/std::pow(sigma_ypHE3Ev3int_G2_SOS, 2);
            break;
        case 391:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE1HE3int_G2_SOS, 2)) + 0.5*std::pow(-mypE1HE3int_G2_SOS + ypE1HE3int_G2_SOS, 2)/std::pow(sigma_ypE1HE3int_G2_SOS, 2);
            break;
        case 392:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE3E4int_G2_SOS, 2)) + 0.5*std::pow(-mypHE3E4int_G2_SOS + ypHE3E4int_G2_SOS, 2)/std::pow(sigma_ypHE3E4int_G2_SOS, 2);
            break;
        case 393:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE3HE4int_G2_SOS, 2)) + 0.5*std::pow(-mypHE3HE4int_G2_SOS + ypHE3HE4int_G2_SOS, 2)/std::pow(sigma_ypHE3HE4int_G2_SOS, 2);
            break;
        case 394:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE2HE4int_G2_SOS, 2)) + 0.5*std::pow(-mypE2HE4int_G2_SOS + ypE2HE4int_G2_SOS, 2)/std::pow(sigma_ypE2HE4int_G2_SOS, 2);
            break;
        case 395:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE4Ev3int_G2_SOS, 2)) + 0.5*std::pow(-mypHE4Ev3int_G2_SOS + ypHE4Ev3int_G2_SOS, 2)/std::pow(sigma_ypHE4Ev3int_G2_SOS, 2);
            break;
        case 396:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE1HE4int_G2_SOS, 2)) + 0.5*std::pow(-mypE1HE4int_G2_SOS + ypE1HE4int_G2_SOS, 2)/std::pow(sigma_ypE1HE4int_G2_SOS, 2);
            break;
        case 397:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE3HE4int_G2_SOS, 2)) + 0.5*std::pow(-mypE3HE4int_G2_SOS + ypE3HE4int_G2_SOS, 2)/std::pow(sigma_ypE3HE4int_G2_SOS, 2);
            break;
        case 398:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE4E4int_G2_SOS, 2)) + 0.5*std::pow(-mypHE4E4int_G2_SOS + ypHE4E4int_G2_SOS, 2)/std::pow(sigma_ypHE4E4int_G2_SOS, 2);
            break;
        case 399:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE4HE4int_G2_SOS, 2)) + 0.5*std::pow(-mypHE4HE4int_G2_SOS + ypHE4HE4int_G2_SOS, 2)/std::pow(sigma_ypHE4HE4int_G2_SOS, 2);
            break;
        case 400:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHGF_Met_Metint_G2_SOS, 2)) + 0.5*std::pow(-mypHGF_Met_Metint_G2_SOS + ypHGF_Met_Metint_G2_SOS, 2)/std::pow(sigma_ypHGF_Met_Metint_G2_SOS, 2);
            break;
        case 401:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHGF_Met_HGF_Metint_G2_SOS, 2)) + 0.5*std::pow(-mypHGF_Met_HGF_Metint_G2_SOS + ypHGF_Met_HGF_Metint_G2_SOS, 2)/std::pow(sigma_ypHGF_Met_HGF_Metint_G2_SOS, 2);
            break;
        case 402:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypPPrPPrint_G2_SOS, 2)) + 0.5*std::pow(-mypPPrPPrint_G2_SOS + ypPPrPPrint_G2_SOS, 2)/std::pow(sigma_ypPPrPPrint_G2_SOS, 2);
            break;
        case 403:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypPPrPrint_G2_SOS, 2)) + 0.5*std::pow(-mypPPrPrint_G2_SOS + ypPPrPrint_G2_SOS, 2)/std::pow(sigma_ypPPrPrint_G2_SOS, 2);
            break;
        case 404:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypFFrFFrint_G2_SOS, 2)) + 0.5*std::pow(-mypFFrFFrint_G2_SOS + ypFFrFFrint_G2_SOS, 2)/std::pow(sigma_ypFFrFFrint_G2_SOS, 2);
            break;
        case 405:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypFFrFrint_G2_SOS, 2)) + 0.5*std::pow(-mypFFrFrint_G2_SOS + ypFFrFrint_G2_SOS, 2)/std::pow(sigma_ypFFrFrint_G2_SOS, 2);
            break;
        case 406:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypIIrIr_int_IRS_G2_SOS, 2)) + 0.5*std::pow(-mypIIrIr_int_IRS_G2_SOS + ypIIrIr_int_IRS_G2_SOS, 2)/std::pow(sigma_ypIIrIr_int_IRS_G2_SOS, 2);
            break;
        case 407:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypINS_Isr_Isr_int_IRS_G2_SOS, 2)) + 0.5*std::pow(-mypINS_Isr_Isr_int_IRS_G2_SOS + ypINS_Isr_Isr_int_IRS_G2_SOS, 2)/std::pow(sigma_ypINS_Isr_Isr_int_IRS_G2_SOS, 2);
            break;
        case 408:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypIIrIrI_int_IRS_G2_SOS, 2)) + 0.5*std::pow(-mypIIrIrI_int_IRS_G2_SOS + ypIIrIrI_int_IRS_G2_SOS, 2)/std::pow(sigma_ypIIrIrI_int_IRS_G2_SOS, 2);
            break;
        case 409:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypINS_Isr_Isr_INS_int_IRS_G2_SOS, 2)) + 0.5*std::pow(-mypINS_Isr_Isr_INS_int_IRS_G2_SOS + ypINS_Isr_Isr_INS_int_IRS_G2_SOS, 2)/std::pow(sigma_ypINS_Isr_Isr_INS_int_IRS_G2_SOS, 2);
            break;
        case 410:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E2_PLCg, 2)) + 0.5*std::pow(-mypEE1E2_PLCg + ypEE1E2_PLCg, 2)/std::pow(sigma_ypEE1E2_PLCg, 2);
            break;
        case 411:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1Ev3_PLCg, 2)) + 0.5*std::pow(-mypEE1Ev3_PLCg + ypEE1Ev3_PLCg, 2)/std::pow(sigma_ypEE1Ev3_PLCg, 2);
            break;
        case 412:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E1_PLCg, 2)) + 0.5*std::pow(-mypEE1E1_PLCg + ypEE1E1_PLCg, 2)/std::pow(sigma_ypEE1E1_PLCg, 2);
            break;
        case 413:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1EE1_PLCg, 2)) + 0.5*std::pow(-mypEE1EE1_PLCg + ypEE1EE1_PLCg, 2)/std::pow(sigma_ypEE1EE1_PLCg, 2);
            break;
        case 414:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E3_PLCg, 2)) + 0.5*std::pow(-mypEE1E3_PLCg + ypEE1E3_PLCg, 2)/std::pow(sigma_ypEE1E3_PLCg, 2);
            break;
        case 415:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1HE3_PLCg, 2)) + 0.5*std::pow(-mypEE1HE3_PLCg + ypEE1HE3_PLCg, 2)/std::pow(sigma_ypEE1HE3_PLCg, 2);
            break;
        case 416:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E4_PLCg, 2)) + 0.5*std::pow(-mypEE1E4_PLCg + ypEE1E4_PLCg, 2)/std::pow(sigma_ypEE1E4_PLCg, 2);
            break;
        case 417:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1HE4_PLCg, 2)) + 0.5*std::pow(-mypEE1HE4_PLCg + ypEE1HE4_PLCg, 2)/std::pow(sigma_ypEE1HE4_PLCg, 2);
            break;
        case 418:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE2HE3_PLCg, 2)) + 0.5*std::pow(-mypE2HE3_PLCg + ypE2HE3_PLCg, 2)/std::pow(sigma_ypE2HE3_PLCg, 2);
            break;
        case 419:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE3Ev3_PLCg, 2)) + 0.5*std::pow(-mypHE3Ev3_PLCg + ypHE3Ev3_PLCg, 2)/std::pow(sigma_ypHE3Ev3_PLCg, 2);
            break;
        case 420:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE1HE3_PLCg, 2)) + 0.5*std::pow(-mypE1HE3_PLCg + ypE1HE3_PLCg, 2)/std::pow(sigma_ypE1HE3_PLCg, 2);
            break;
        case 421:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE3E4_PLCg, 2)) + 0.5*std::pow(-mypHE3E4_PLCg + ypHE3E4_PLCg, 2)/std::pow(sigma_ypHE3E4_PLCg, 2);
            break;
        case 422:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE3HE4_PLCg, 2)) + 0.5*std::pow(-mypHE3HE4_PLCg + ypHE3HE4_PLCg, 2)/std::pow(sigma_ypHE3HE4_PLCg, 2);
            break;
        case 423:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE2HE4_PLCg, 2)) + 0.5*std::pow(-mypE2HE4_PLCg + ypE2HE4_PLCg, 2)/std::pow(sigma_ypE2HE4_PLCg, 2);
            break;
        case 424:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE4Ev3_PLCg, 2)) + 0.5*std::pow(-mypHE4Ev3_PLCg + ypHE4Ev3_PLCg, 2)/std::pow(sigma_ypHE4Ev3_PLCg, 2);
            break;
        case 425:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE1HE4_PLCg, 2)) + 0.5*std::pow(-mypE1HE4_PLCg + ypE1HE4_PLCg, 2)/std::pow(sigma_ypE1HE4_PLCg, 2);
            break;
        case 426:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE3HE4_PLCg, 2)) + 0.5*std::pow(-mypE3HE4_PLCg + ypE3HE4_PLCg, 2)/std::pow(sigma_ypE3HE4_PLCg, 2);
            break;
        case 427:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE4E4_PLCg, 2)) + 0.5*std::pow(-mypHE4E4_PLCg + ypHE4E4_PLCg, 2)/std::pow(sigma_ypHE4E4_PLCg, 2);
            break;
        case 428:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE4HE4_PLCg, 2)) + 0.5*std::pow(-mypHE4HE4_PLCg + ypHE4HE4_PLCg, 2)/std::pow(sigma_ypHE4HE4_PLCg, 2);
            break;
        case 429:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHGF_Met_Met_PLCg, 2)) + 0.5*std::pow(-mypHGF_Met_Met_PLCg + ypHGF_Met_Met_PLCg, 2)/std::pow(sigma_ypHGF_Met_Met_PLCg, 2);
            break;
        case 430:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHGF_Met_HGF_Met_PLCg, 2)) + 0.5*std::pow(-mypHGF_Met_HGF_Met_PLCg + ypHGF_Met_HGF_Met_PLCg, 2)/std::pow(sigma_ypHGF_Met_HGF_Met_PLCg, 2);
            break;
        case 431:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypPPrPPr_PLCg, 2)) + 0.5*std::pow(-mypPPrPPr_PLCg + ypPPrPPr_PLCg, 2)/std::pow(sigma_ypPPrPPr_PLCg, 2);
            break;
        case 432:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypPPrPr_PLCg, 2)) + 0.5*std::pow(-mypPPrPr_PLCg + ypPPrPr_PLCg, 2)/std::pow(sigma_ypPPrPr_PLCg, 2);
            break;
        case 433:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypFFrFFr_PLCg, 2)) + 0.5*std::pow(-mypFFrFFr_PLCg + ypFFrFFr_PLCg, 2)/std::pow(sigma_ypFFrFFr_PLCg, 2);
            break;
        case 434:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypFFrFr_PLCg, 2)) + 0.5*std::pow(-mypFFrFr_PLCg + ypFFrFr_PLCg, 2)/std::pow(sigma_ypFFrFr_PLCg, 2);
            break;
        case 435:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypIIrIr_IRS_PLCg, 2)) + 0.5*std::pow(-mypIIrIr_IRS_PLCg + ypIIrIr_IRS_PLCg, 2)/std::pow(sigma_ypIIrIr_IRS_PLCg, 2);
            break;
        case 436:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypINS_Isr_Isr_IRS_PLCg, 2)) + 0.5*std::pow(-mypINS_Isr_Isr_IRS_PLCg + ypINS_Isr_Isr_IRS_PLCg, 2)/std::pow(sigma_ypINS_Isr_Isr_IRS_PLCg, 2);
            break;
        case 437:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypIIrIrI_IRS_PLCg, 2)) + 0.5*std::pow(-mypIIrIrI_IRS_PLCg + ypIIrIrI_IRS_PLCg, 2)/std::pow(sigma_ypIIrIrI_IRS_PLCg, 2);
            break;
        case 438:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypINS_Isr_Isr_INS_IRS_PLCg, 2)) + 0.5*std::pow(-mypINS_Isr_Isr_INS_IRS_PLCg + ypINS_Isr_Isr_INS_IRS_PLCg, 2)/std::pow(sigma_ypINS_Isr_Isr_INS_IRS_PLCg, 2);
            break;
        case 439:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E2_PI3K1, 2)) + 0.5*std::pow(-mypEE1E2_PI3K1 + ypEE1E2_PI3K1, 2)/std::pow(sigma_ypEE1E2_PI3K1, 2);
            break;
        case 440:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1Ev3_PI3K1, 2)) + 0.5*std::pow(-mypEE1Ev3_PI3K1 + ypEE1Ev3_PI3K1, 2)/std::pow(sigma_ypEE1Ev3_PI3K1, 2);
            break;
        case 441:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E1_PI3K1, 2)) + 0.5*std::pow(-mypEE1E1_PI3K1 + ypEE1E1_PI3K1, 2)/std::pow(sigma_ypEE1E1_PI3K1, 2);
            break;
        case 442:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1EE1_PI3K1, 2)) + 0.5*std::pow(-mypEE1EE1_PI3K1 + ypEE1EE1_PI3K1, 2)/std::pow(sigma_ypEE1EE1_PI3K1, 2);
            break;
        case 443:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E3_PI3K1, 2)) + 0.5*std::pow(-mypEE1E3_PI3K1 + ypEE1E3_PI3K1, 2)/std::pow(sigma_ypEE1E3_PI3K1, 2);
            break;
        case 444:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1HE3_PI3K1, 2)) + 0.5*std::pow(-mypEE1HE3_PI3K1 + ypEE1HE3_PI3K1, 2)/std::pow(sigma_ypEE1HE3_PI3K1, 2);
            break;
        case 445:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E4_PI3K1, 2)) + 0.5*std::pow(-mypEE1E4_PI3K1 + ypEE1E4_PI3K1, 2)/std::pow(sigma_ypEE1E4_PI3K1, 2);
            break;
        case 446:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1HE4_PI3K1, 2)) + 0.5*std::pow(-mypEE1HE4_PI3K1 + ypEE1HE4_PI3K1, 2)/std::pow(sigma_ypEE1HE4_PI3K1, 2);
            break;
        case 447:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE2HE3_PI3K1, 2)) + 0.5*std::pow(-mypE2HE3_PI3K1 + ypE2HE3_PI3K1, 2)/std::pow(sigma_ypE2HE3_PI3K1, 2);
            break;
        case 448:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE3Ev3_PI3K1, 2)) + 0.5*std::pow(-mypHE3Ev3_PI3K1 + ypHE3Ev3_PI3K1, 2)/std::pow(sigma_ypHE3Ev3_PI3K1, 2);
            break;
        case 449:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE1HE3_PI3K1, 2)) + 0.5*std::pow(-mypE1HE3_PI3K1 + ypE1HE3_PI3K1, 2)/std::pow(sigma_ypE1HE3_PI3K1, 2);
            break;
        case 450:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE3E4_PI3K1, 2)) + 0.5*std::pow(-mypHE3E4_PI3K1 + ypHE3E4_PI3K1, 2)/std::pow(sigma_ypHE3E4_PI3K1, 2);
            break;
        case 451:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE3HE4_PI3K1, 2)) + 0.5*std::pow(-mypHE3HE4_PI3K1 + ypHE3HE4_PI3K1, 2)/std::pow(sigma_ypHE3HE4_PI3K1, 2);
            break;
        case 452:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE2HE4_PI3K1, 2)) + 0.5*std::pow(-mypE2HE4_PI3K1 + ypE2HE4_PI3K1, 2)/std::pow(sigma_ypE2HE4_PI3K1, 2);
            break;
        case 453:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE4Ev3_PI3K1, 2)) + 0.5*std::pow(-mypHE4Ev3_PI3K1 + ypHE4Ev3_PI3K1, 2)/std::pow(sigma_ypHE4Ev3_PI3K1, 2);
            break;
        case 454:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE1HE4_PI3K1, 2)) + 0.5*std::pow(-mypE1HE4_PI3K1 + ypE1HE4_PI3K1, 2)/std::pow(sigma_ypE1HE4_PI3K1, 2);
            break;
        case 455:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE3HE4_PI3K1, 2)) + 0.5*std::pow(-mypE3HE4_PI3K1 + ypE3HE4_PI3K1, 2)/std::pow(sigma_ypE3HE4_PI3K1, 2);
            break;
        case 456:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE4E4_PI3K1, 2)) + 0.5*std::pow(-mypHE4E4_PI3K1 + ypHE4E4_PI3K1, 2)/std::pow(sigma_ypHE4E4_PI3K1, 2);
            break;
        case 457:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE4HE4_PI3K1, 2)) + 0.5*std::pow(-mypHE4HE4_PI3K1 + ypHE4HE4_PI3K1, 2)/std::pow(sigma_ypHE4HE4_PI3K1, 2);
            break;
        case 458:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHGF_Met_Met_PI3K1, 2)) + 0.5*std::pow(-mypHGF_Met_Met_PI3K1 + ypHGF_Met_Met_PI3K1, 2)/std::pow(sigma_ypHGF_Met_Met_PI3K1, 2);
            break;
        case 459:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHGF_Met_HGF_Met_PI3K1, 2)) + 0.5*std::pow(-mypHGF_Met_HGF_Met_PI3K1 + ypHGF_Met_HGF_Met_PI3K1, 2)/std::pow(sigma_ypHGF_Met_HGF_Met_PI3K1, 2);
            break;
        case 460:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypPPrPPr_PI3K1, 2)) + 0.5*std::pow(-mypPPrPPr_PI3K1 + ypPPrPPr_PI3K1, 2)/std::pow(sigma_ypPPrPPr_PI3K1, 2);
            break;
        case 461:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypPPrPr_PI3K1, 2)) + 0.5*std::pow(-mypPPrPr_PI3K1 + ypPPrPr_PI3K1, 2)/std::pow(sigma_ypPPrPr_PI3K1, 2);
            break;
        case 462:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypFFrFFr_PI3K1, 2)) + 0.5*std::pow(-mypFFrFFr_PI3K1 + ypFFrFFr_PI3K1, 2)/std::pow(sigma_ypFFrFFr_PI3K1, 2);
            break;
        case 463:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypFFrFr_PI3K1, 2)) + 0.5*std::pow(-mypFFrFr_PI3K1 + ypFFrFr_PI3K1, 2)/std::pow(sigma_ypFFrFr_PI3K1, 2);
            break;
        case 464:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypIIrIr_IRS_PI3K1, 2)) + 0.5*std::pow(-mypIIrIr_IRS_PI3K1 + ypIIrIr_IRS_PI3K1, 2)/std::pow(sigma_ypIIrIr_IRS_PI3K1, 2);
            break;
        case 465:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypINS_Isr_Isr_IRS_PI3K1, 2)) + 0.5*std::pow(-mypINS_Isr_Isr_IRS_PI3K1 + ypINS_Isr_Isr_IRS_PI3K1, 2)/std::pow(sigma_ypINS_Isr_Isr_IRS_PI3K1, 2);
            break;
        case 466:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypIIrIrI_IRS_PI3K1, 2)) + 0.5*std::pow(-mypIIrIrI_IRS_PI3K1 + ypIIrIrI_IRS_PI3K1, 2)/std::pow(sigma_ypIIrIrI_IRS_PI3K1, 2);
            break;
        case 467:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypINS_Isr_Isr_INS_IRS_PI3K1, 2)) + 0.5*std::pow(-mypINS_Isr_Isr_INS_IRS_PI3K1 + ypINS_Isr_Isr_INS_IRS_PI3K1, 2)/std::pow(sigma_ypINS_Isr_Isr_INS_IRS_PI3K1, 2);
            break;
        case 468:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E2_PI3K2, 2)) + 0.5*std::pow(-mypEE1E2_PI3K2 + ypEE1E2_PI3K2, 2)/std::pow(sigma_ypEE1E2_PI3K2, 2);
            break;
        case 469:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1Ev3_PI3K2, 2)) + 0.5*std::pow(-mypEE1Ev3_PI3K2 + ypEE1Ev3_PI3K2, 2)/std::pow(sigma_ypEE1Ev3_PI3K2, 2);
            break;
        case 470:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E1_PI3K2, 2)) + 0.5*std::pow(-mypEE1E1_PI3K2 + ypEE1E1_PI3K2, 2)/std::pow(sigma_ypEE1E1_PI3K2, 2);
            break;
        case 471:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1EE1_PI3K2, 2)) + 0.5*std::pow(-mypEE1EE1_PI3K2 + ypEE1EE1_PI3K2, 2)/std::pow(sigma_ypEE1EE1_PI3K2, 2);
            break;
        case 472:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E3_PI3K2, 2)) + 0.5*std::pow(-mypEE1E3_PI3K2 + ypEE1E3_PI3K2, 2)/std::pow(sigma_ypEE1E3_PI3K2, 2);
            break;
        case 473:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1HE3_PI3K2, 2)) + 0.5*std::pow(-mypEE1HE3_PI3K2 + ypEE1HE3_PI3K2, 2)/std::pow(sigma_ypEE1HE3_PI3K2, 2);
            break;
        case 474:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E4_PI3K2, 2)) + 0.5*std::pow(-mypEE1E4_PI3K2 + ypEE1E4_PI3K2, 2)/std::pow(sigma_ypEE1E4_PI3K2, 2);
            break;
        case 475:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1HE4_PI3K2, 2)) + 0.5*std::pow(-mypEE1HE4_PI3K2 + ypEE1HE4_PI3K2, 2)/std::pow(sigma_ypEE1HE4_PI3K2, 2);
            break;
        case 476:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE2HE3_PI3K2, 2)) + 0.5*std::pow(-mypE2HE3_PI3K2 + ypE2HE3_PI3K2, 2)/std::pow(sigma_ypE2HE3_PI3K2, 2);
            break;
        case 477:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE3Ev3_PI3K2, 2)) + 0.5*std::pow(-mypHE3Ev3_PI3K2 + ypHE3Ev3_PI3K2, 2)/std::pow(sigma_ypHE3Ev3_PI3K2, 2);
            break;
        case 478:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE1HE3_PI3K2, 2)) + 0.5*std::pow(-mypE1HE3_PI3K2 + ypE1HE3_PI3K2, 2)/std::pow(sigma_ypE1HE3_PI3K2, 2);
            break;
        case 479:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE3E4_PI3K2, 2)) + 0.5*std::pow(-mypHE3E4_PI3K2 + ypHE3E4_PI3K2, 2)/std::pow(sigma_ypHE3E4_PI3K2, 2);
            break;
        case 480:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE3HE4_PI3K2, 2)) + 0.5*std::pow(-mypHE3HE4_PI3K2 + ypHE3HE4_PI3K2, 2)/std::pow(sigma_ypHE3HE4_PI3K2, 2);
            break;
        case 481:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE2HE4_PI3K2, 2)) + 0.5*std::pow(-mypE2HE4_PI3K2 + ypE2HE4_PI3K2, 2)/std::pow(sigma_ypE2HE4_PI3K2, 2);
            break;
        case 482:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE4Ev3_PI3K2, 2)) + 0.5*std::pow(-mypHE4Ev3_PI3K2 + ypHE4Ev3_PI3K2, 2)/std::pow(sigma_ypHE4Ev3_PI3K2, 2);
            break;
        case 483:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE1HE4_PI3K2, 2)) + 0.5*std::pow(-mypE1HE4_PI3K2 + ypE1HE4_PI3K2, 2)/std::pow(sigma_ypE1HE4_PI3K2, 2);
            break;
        case 484:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE3HE4_PI3K2, 2)) + 0.5*std::pow(-mypE3HE4_PI3K2 + ypE3HE4_PI3K2, 2)/std::pow(sigma_ypE3HE4_PI3K2, 2);
            break;
        case 485:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE4E4_PI3K2, 2)) + 0.5*std::pow(-mypHE4E4_PI3K2 + ypHE4E4_PI3K2, 2)/std::pow(sigma_ypHE4E4_PI3K2, 2);
            break;
        case 486:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE4HE4_PI3K2, 2)) + 0.5*std::pow(-mypHE4HE4_PI3K2 + ypHE4HE4_PI3K2, 2)/std::pow(sigma_ypHE4HE4_PI3K2, 2);
            break;
        case 487:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHGF_Met_Met_PI3K2, 2)) + 0.5*std::pow(-mypHGF_Met_Met_PI3K2 + ypHGF_Met_Met_PI3K2, 2)/std::pow(sigma_ypHGF_Met_Met_PI3K2, 2);
            break;
        case 488:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHGF_Met_HGF_Met_PI3K2, 2)) + 0.5*std::pow(-mypHGF_Met_HGF_Met_PI3K2 + ypHGF_Met_HGF_Met_PI3K2, 2)/std::pow(sigma_ypHGF_Met_HGF_Met_PI3K2, 2);
            break;
        case 489:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypPPrPPr_PI3K2, 2)) + 0.5*std::pow(-mypPPrPPr_PI3K2 + ypPPrPPr_PI3K2, 2)/std::pow(sigma_ypPPrPPr_PI3K2, 2);
            break;
        case 490:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypPPrPr_PI3K2, 2)) + 0.5*std::pow(-mypPPrPr_PI3K2 + ypPPrPr_PI3K2, 2)/std::pow(sigma_ypPPrPr_PI3K2, 2);
            break;
        case 491:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypFFrFFr_PI3K2, 2)) + 0.5*std::pow(-mypFFrFFr_PI3K2 + ypFFrFFr_PI3K2, 2)/std::pow(sigma_ypFFrFFr_PI3K2, 2);
            break;
        case 492:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypFFrFr_PI3K2, 2)) + 0.5*std::pow(-mypFFrFr_PI3K2 + ypFFrFr_PI3K2, 2)/std::pow(sigma_ypFFrFr_PI3K2, 2);
            break;
        case 493:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypIIrIr_IRS_PI3K2, 2)) + 0.5*std::pow(-mypIIrIr_IRS_PI3K2 + ypIIrIr_IRS_PI3K2, 2)/std::pow(sigma_ypIIrIr_IRS_PI3K2, 2);
            break;
        case 494:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypINS_Isr_Isr_IRS_PI3K2, 2)) + 0.5*std::pow(-mypINS_Isr_Isr_IRS_PI3K2 + ypINS_Isr_Isr_IRS_PI3K2, 2)/std::pow(sigma_ypINS_Isr_Isr_IRS_PI3K2, 2);
            break;
        case 495:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypIIrIrI_IRS_PI3K2, 2)) + 0.5*std::pow(-mypIIrIrI_IRS_PI3K2 + ypIIrIrI_IRS_PI3K2, 2)/std::pow(sigma_ypIIrIrI_IRS_PI3K2, 2);
            break;
        case 496:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypINS_Isr_Isr_INS_IRS_PI3K2, 2)) + 0.5*std::pow(-mypINS_Isr_Isr_INS_IRS_PI3K2 + ypINS_Isr_Isr_INS_IRS_PI3K2, 2)/std::pow(sigma_ypINS_Isr_Isr_INS_IRS_PI3K2, 2);
            break;
        case 497:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E2int_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypEE1E2int_G2_SOS_RasD + ypEE1E2int_G2_SOS_RasD, 2)/std::pow(sigma_ypEE1E2int_G2_SOS_RasD, 2);
            break;
        case 498:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1Ev3int_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypEE1Ev3int_G2_SOS_RasD + ypEE1Ev3int_G2_SOS_RasD, 2)/std::pow(sigma_ypEE1Ev3int_G2_SOS_RasD, 2);
            break;
        case 499:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E1int_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypEE1E1int_G2_SOS_RasD + ypEE1E1int_G2_SOS_RasD, 2)/std::pow(sigma_ypEE1E1int_G2_SOS_RasD, 2);
            break;
        case 500:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1EE1int_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypEE1EE1int_G2_SOS_RasD + ypEE1EE1int_G2_SOS_RasD, 2)/std::pow(sigma_ypEE1EE1int_G2_SOS_RasD, 2);
            break;
        case 501:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E3int_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypEE1E3int_G2_SOS_RasD + ypEE1E3int_G2_SOS_RasD, 2)/std::pow(sigma_ypEE1E3int_G2_SOS_RasD, 2);
            break;
        case 502:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1HE3int_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypEE1HE3int_G2_SOS_RasD + ypEE1HE3int_G2_SOS_RasD, 2)/std::pow(sigma_ypEE1HE3int_G2_SOS_RasD, 2);
            break;
        case 503:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E4int_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypEE1E4int_G2_SOS_RasD + ypEE1E4int_G2_SOS_RasD, 2)/std::pow(sigma_ypEE1E4int_G2_SOS_RasD, 2);
            break;
        case 504:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1HE4int_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypEE1HE4int_G2_SOS_RasD + ypEE1HE4int_G2_SOS_RasD, 2)/std::pow(sigma_ypEE1HE4int_G2_SOS_RasD, 2);
            break;
        case 505:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE2HE3int_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypE2HE3int_G2_SOS_RasD + ypE2HE3int_G2_SOS_RasD, 2)/std::pow(sigma_ypE2HE3int_G2_SOS_RasD, 2);
            break;
        case 506:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE3Ev3int_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypHE3Ev3int_G2_SOS_RasD + ypHE3Ev3int_G2_SOS_RasD, 2)/std::pow(sigma_ypHE3Ev3int_G2_SOS_RasD, 2);
            break;
        case 507:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE1HE3int_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypE1HE3int_G2_SOS_RasD + ypE1HE3int_G2_SOS_RasD, 2)/std::pow(sigma_ypE1HE3int_G2_SOS_RasD, 2);
            break;
        case 508:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE3E4int_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypHE3E4int_G2_SOS_RasD + ypHE3E4int_G2_SOS_RasD, 2)/std::pow(sigma_ypHE3E4int_G2_SOS_RasD, 2);
            break;
        case 509:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE3HE4int_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypHE3HE4int_G2_SOS_RasD + ypHE3HE4int_G2_SOS_RasD, 2)/std::pow(sigma_ypHE3HE4int_G2_SOS_RasD, 2);
            break;
        case 510:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE2HE4int_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypE2HE4int_G2_SOS_RasD + ypE2HE4int_G2_SOS_RasD, 2)/std::pow(sigma_ypE2HE4int_G2_SOS_RasD, 2);
            break;
        case 511:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE4Ev3int_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypHE4Ev3int_G2_SOS_RasD + ypHE4Ev3int_G2_SOS_RasD, 2)/std::pow(sigma_ypHE4Ev3int_G2_SOS_RasD, 2);
            break;
        case 512:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE1HE4int_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypE1HE4int_G2_SOS_RasD + ypE1HE4int_G2_SOS_RasD, 2)/std::pow(sigma_ypE1HE4int_G2_SOS_RasD, 2);
            break;
        case 513:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE3HE4int_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypE3HE4int_G2_SOS_RasD + ypE3HE4int_G2_SOS_RasD, 2)/std::pow(sigma_ypE3HE4int_G2_SOS_RasD, 2);
            break;
        case 514:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE4E4int_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypHE4E4int_G2_SOS_RasD + ypHE4E4int_G2_SOS_RasD, 2)/std::pow(sigma_ypHE4E4int_G2_SOS_RasD, 2);
            break;
        case 515:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE4HE4int_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypHE4HE4int_G2_SOS_RasD + ypHE4HE4int_G2_SOS_RasD, 2)/std::pow(sigma_ypHE4HE4int_G2_SOS_RasD, 2);
            break;
        case 516:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHGF_Met_Metint_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypHGF_Met_Metint_G2_SOS_RasD + ypHGF_Met_Metint_G2_SOS_RasD, 2)/std::pow(sigma_ypHGF_Met_Metint_G2_SOS_RasD, 2);
            break;
        case 517:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHGF_Met_HGF_Metint_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypHGF_Met_HGF_Metint_G2_SOS_RasD + ypHGF_Met_HGF_Metint_G2_SOS_RasD, 2)/std::pow(sigma_ypHGF_Met_HGF_Metint_G2_SOS_RasD, 2);
            break;
        case 518:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypPPrPPrint_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypPPrPPrint_G2_SOS_RasD + ypPPrPPrint_G2_SOS_RasD, 2)/std::pow(sigma_ypPPrPPrint_G2_SOS_RasD, 2);
            break;
        case 519:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypPPrPrint_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypPPrPrint_G2_SOS_RasD + ypPPrPrint_G2_SOS_RasD, 2)/std::pow(sigma_ypPPrPrint_G2_SOS_RasD, 2);
            break;
        case 520:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypFFrFFrint_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypFFrFFrint_G2_SOS_RasD + ypFFrFFrint_G2_SOS_RasD, 2)/std::pow(sigma_ypFFrFFrint_G2_SOS_RasD, 2);
            break;
        case 521:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypFFrFrint_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypFFrFrint_G2_SOS_RasD + ypFFrFrint_G2_SOS_RasD, 2)/std::pow(sigma_ypFFrFrint_G2_SOS_RasD, 2);
            break;
        case 522:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypIIrIr_int_IRS_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypIIrIr_int_IRS_G2_SOS_RasD + ypIIrIr_int_IRS_G2_SOS_RasD, 2)/std::pow(sigma_ypIIrIr_int_IRS_G2_SOS_RasD, 2);
            break;
        case 523:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypINS_Isr_Isr_int_IRS_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypINS_Isr_Isr_int_IRS_G2_SOS_RasD + ypINS_Isr_Isr_int_IRS_G2_SOS_RasD, 2)/std::pow(sigma_ypINS_Isr_Isr_int_IRS_G2_SOS_RasD, 2);
            break;
        case 524:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypIIrIrI_int_IRS_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypIIrIrI_int_IRS_G2_SOS_RasD + ypIIrIrI_int_IRS_G2_SOS_RasD, 2)/std::pow(sigma_ypIIrIrI_int_IRS_G2_SOS_RasD, 2);
            break;
        case 525:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypINS_Isr_Isr_INS_int_IRS_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypINS_Isr_Isr_INS_int_IRS_G2_SOS_RasD + ypINS_Isr_Isr_INS_int_IRS_G2_SOS_RasD, 2)/std::pow(sigma_ypINS_Isr_Isr_INS_int_IRS_G2_SOS_RasD, 2);
            break;
        case 526:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E2_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypEE1E2_G2_SOS_RasD + ypEE1E2_G2_SOS_RasD, 2)/std::pow(sigma_ypEE1E2_G2_SOS_RasD, 2);
            break;
        case 527:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1Ev3_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypEE1Ev3_G2_SOS_RasD + ypEE1Ev3_G2_SOS_RasD, 2)/std::pow(sigma_ypEE1Ev3_G2_SOS_RasD, 2);
            break;
        case 528:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E1_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypEE1E1_G2_SOS_RasD + ypEE1E1_G2_SOS_RasD, 2)/std::pow(sigma_ypEE1E1_G2_SOS_RasD, 2);
            break;
        case 529:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1EE1_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypEE1EE1_G2_SOS_RasD + ypEE1EE1_G2_SOS_RasD, 2)/std::pow(sigma_ypEE1EE1_G2_SOS_RasD, 2);
            break;
        case 530:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E3_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypEE1E3_G2_SOS_RasD + ypEE1E3_G2_SOS_RasD, 2)/std::pow(sigma_ypEE1E3_G2_SOS_RasD, 2);
            break;
        case 531:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1HE3_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypEE1HE3_G2_SOS_RasD + ypEE1HE3_G2_SOS_RasD, 2)/std::pow(sigma_ypEE1HE3_G2_SOS_RasD, 2);
            break;
        case 532:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E4_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypEE1E4_G2_SOS_RasD + ypEE1E4_G2_SOS_RasD, 2)/std::pow(sigma_ypEE1E4_G2_SOS_RasD, 2);
            break;
        case 533:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1HE4_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypEE1HE4_G2_SOS_RasD + ypEE1HE4_G2_SOS_RasD, 2)/std::pow(sigma_ypEE1HE4_G2_SOS_RasD, 2);
            break;
        case 534:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE2HE3_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypE2HE3_G2_SOS_RasD + ypE2HE3_G2_SOS_RasD, 2)/std::pow(sigma_ypE2HE3_G2_SOS_RasD, 2);
            break;
        case 535:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE3Ev3_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypHE3Ev3_G2_SOS_RasD + ypHE3Ev3_G2_SOS_RasD, 2)/std::pow(sigma_ypHE3Ev3_G2_SOS_RasD, 2);
            break;
        case 536:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE1HE3_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypE1HE3_G2_SOS_RasD + ypE1HE3_G2_SOS_RasD, 2)/std::pow(sigma_ypE1HE3_G2_SOS_RasD, 2);
            break;
        case 537:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE3E4_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypHE3E4_G2_SOS_RasD + ypHE3E4_G2_SOS_RasD, 2)/std::pow(sigma_ypHE3E4_G2_SOS_RasD, 2);
            break;
        case 538:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE3HE4_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypHE3HE4_G2_SOS_RasD + ypHE3HE4_G2_SOS_RasD, 2)/std::pow(sigma_ypHE3HE4_G2_SOS_RasD, 2);
            break;
        case 539:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE2HE4_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypE2HE4_G2_SOS_RasD + ypE2HE4_G2_SOS_RasD, 2)/std::pow(sigma_ypE2HE4_G2_SOS_RasD, 2);
            break;
        case 540:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE4Ev3_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypHE4Ev3_G2_SOS_RasD + ypHE4Ev3_G2_SOS_RasD, 2)/std::pow(sigma_ypHE4Ev3_G2_SOS_RasD, 2);
            break;
        case 541:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE1HE4_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypE1HE4_G2_SOS_RasD + ypE1HE4_G2_SOS_RasD, 2)/std::pow(sigma_ypE1HE4_G2_SOS_RasD, 2);
            break;
        case 542:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE3HE4_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypE3HE4_G2_SOS_RasD + ypE3HE4_G2_SOS_RasD, 2)/std::pow(sigma_ypE3HE4_G2_SOS_RasD, 2);
            break;
        case 543:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE4E4_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypHE4E4_G2_SOS_RasD + ypHE4E4_G2_SOS_RasD, 2)/std::pow(sigma_ypHE4E4_G2_SOS_RasD, 2);
            break;
        case 544:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE4HE4_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypHE4HE4_G2_SOS_RasD + ypHE4HE4_G2_SOS_RasD, 2)/std::pow(sigma_ypHE4HE4_G2_SOS_RasD, 2);
            break;
        case 545:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHGF_Met_Met_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypHGF_Met_Met_G2_SOS_RasD + ypHGF_Met_Met_G2_SOS_RasD, 2)/std::pow(sigma_ypHGF_Met_Met_G2_SOS_RasD, 2);
            break;
        case 546:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHGF_Met_HGF_Met_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypHGF_Met_HGF_Met_G2_SOS_RasD + ypHGF_Met_HGF_Met_G2_SOS_RasD, 2)/std::pow(sigma_ypHGF_Met_HGF_Met_G2_SOS_RasD, 2);
            break;
        case 547:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypPPrPPr_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypPPrPPr_G2_SOS_RasD + ypPPrPPr_G2_SOS_RasD, 2)/std::pow(sigma_ypPPrPPr_G2_SOS_RasD, 2);
            break;
        case 548:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypPPrPr_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypPPrPr_G2_SOS_RasD + ypPPrPr_G2_SOS_RasD, 2)/std::pow(sigma_ypPPrPr_G2_SOS_RasD, 2);
            break;
        case 549:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypFFrFFr_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypFFrFFr_G2_SOS_RasD + ypFFrFFr_G2_SOS_RasD, 2)/std::pow(sigma_ypFFrFFr_G2_SOS_RasD, 2);
            break;
        case 550:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypFFrFr_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypFFrFr_G2_SOS_RasD + ypFFrFr_G2_SOS_RasD, 2)/std::pow(sigma_ypFFrFr_G2_SOS_RasD, 2);
            break;
        case 551:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypIIrIr_IRS_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypIIrIr_IRS_G2_SOS_RasD + ypIIrIr_IRS_G2_SOS_RasD, 2)/std::pow(sigma_ypIIrIr_IRS_G2_SOS_RasD, 2);
            break;
        case 552:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypINS_Isr_Isr_IRS_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypINS_Isr_Isr_IRS_G2_SOS_RasD + ypINS_Isr_Isr_IRS_G2_SOS_RasD, 2)/std::pow(sigma_ypINS_Isr_Isr_IRS_G2_SOS_RasD, 2);
            break;
        case 553:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypIIrIrI_IRS_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypIIrIrI_IRS_G2_SOS_RasD + ypIIrIrI_IRS_G2_SOS_RasD, 2)/std::pow(sigma_ypIIrIrI_IRS_G2_SOS_RasD, 2);
            break;
        case 554:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypINS_Isr_Isr_INS_IRS_G2_SOS_RasD, 2)) + 0.5*std::pow(-mypINS_Isr_Isr_INS_IRS_G2_SOS_RasD + ypINS_Isr_Isr_INS_IRS_G2_SOS_RasD, 2)/std::pow(sigma_ypINS_Isr_Isr_INS_IRS_G2_SOS_RasD, 2);
            break;
        case 555:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E2_PLCg_PIP2, 2)) + 0.5*std::pow(-mypEE1E2_PLCg_PIP2 + ypEE1E2_PLCg_PIP2, 2)/std::pow(sigma_ypEE1E2_PLCg_PIP2, 2);
            break;
        case 556:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1Ev3_PLCg_PIP2, 2)) + 0.5*std::pow(-mypEE1Ev3_PLCg_PIP2 + ypEE1Ev3_PLCg_PIP2, 2)/std::pow(sigma_ypEE1Ev3_PLCg_PIP2, 2);
            break;
        case 557:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E1_PLCg_PIP2, 2)) + 0.5*std::pow(-mypEE1E1_PLCg_PIP2 + ypEE1E1_PLCg_PIP2, 2)/std::pow(sigma_ypEE1E1_PLCg_PIP2, 2);
            break;
        case 558:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1EE1_PLCg_PIP2, 2)) + 0.5*std::pow(-mypEE1EE1_PLCg_PIP2 + ypEE1EE1_PLCg_PIP2, 2)/std::pow(sigma_ypEE1EE1_PLCg_PIP2, 2);
            break;
        case 559:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E3_PLCg_PIP2, 2)) + 0.5*std::pow(-mypEE1E3_PLCg_PIP2 + ypEE1E3_PLCg_PIP2, 2)/std::pow(sigma_ypEE1E3_PLCg_PIP2, 2);
            break;
        case 560:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1HE3_PLCg_PIP2, 2)) + 0.5*std::pow(-mypEE1HE3_PLCg_PIP2 + ypEE1HE3_PLCg_PIP2, 2)/std::pow(sigma_ypEE1HE3_PLCg_PIP2, 2);
            break;
        case 561:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E4_PLCg_PIP2, 2)) + 0.5*std::pow(-mypEE1E4_PLCg_PIP2 + ypEE1E4_PLCg_PIP2, 2)/std::pow(sigma_ypEE1E4_PLCg_PIP2, 2);
            break;
        case 562:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1HE4_PLCg_PIP2, 2)) + 0.5*std::pow(-mypEE1HE4_PLCg_PIP2 + ypEE1HE4_PLCg_PIP2, 2)/std::pow(sigma_ypEE1HE4_PLCg_PIP2, 2);
            break;
        case 563:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE2HE3_PLCg_PIP2, 2)) + 0.5*std::pow(-mypE2HE3_PLCg_PIP2 + ypE2HE3_PLCg_PIP2, 2)/std::pow(sigma_ypE2HE3_PLCg_PIP2, 2);
            break;
        case 564:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE3Ev3_PLCg_PIP2, 2)) + 0.5*std::pow(-mypHE3Ev3_PLCg_PIP2 + ypHE3Ev3_PLCg_PIP2, 2)/std::pow(sigma_ypHE3Ev3_PLCg_PIP2, 2);
            break;
        case 565:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE1HE3_PLCg_PIP2, 2)) + 0.5*std::pow(-mypE1HE3_PLCg_PIP2 + ypE1HE3_PLCg_PIP2, 2)/std::pow(sigma_ypE1HE3_PLCg_PIP2, 2);
            break;
        case 566:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE3E4_PLCg_PIP2, 2)) + 0.5*std::pow(-mypHE3E4_PLCg_PIP2 + ypHE3E4_PLCg_PIP2, 2)/std::pow(sigma_ypHE3E4_PLCg_PIP2, 2);
            break;
        case 567:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE3HE4_PLCg_PIP2, 2)) + 0.5*std::pow(-mypHE3HE4_PLCg_PIP2 + ypHE3HE4_PLCg_PIP2, 2)/std::pow(sigma_ypHE3HE4_PLCg_PIP2, 2);
            break;
        case 568:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE2HE4_PLCg_PIP2, 2)) + 0.5*std::pow(-mypE2HE4_PLCg_PIP2 + ypE2HE4_PLCg_PIP2, 2)/std::pow(sigma_ypE2HE4_PLCg_PIP2, 2);
            break;
        case 569:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE4Ev3_PLCg_PIP2, 2)) + 0.5*std::pow(-mypHE4Ev3_PLCg_PIP2 + ypHE4Ev3_PLCg_PIP2, 2)/std::pow(sigma_ypHE4Ev3_PLCg_PIP2, 2);
            break;
        case 570:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE1HE4_PLCg_PIP2, 2)) + 0.5*std::pow(-mypE1HE4_PLCg_PIP2 + ypE1HE4_PLCg_PIP2, 2)/std::pow(sigma_ypE1HE4_PLCg_PIP2, 2);
            break;
        case 571:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE3HE4_PLCg_PIP2, 2)) + 0.5*std::pow(-mypE3HE4_PLCg_PIP2 + ypE3HE4_PLCg_PIP2, 2)/std::pow(sigma_ypE3HE4_PLCg_PIP2, 2);
            break;
        case 572:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE4E4_PLCg_PIP2, 2)) + 0.5*std::pow(-mypHE4E4_PLCg_PIP2 + ypHE4E4_PLCg_PIP2, 2)/std::pow(sigma_ypHE4E4_PLCg_PIP2, 2);
            break;
        case 573:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE4HE4_PLCg_PIP2, 2)) + 0.5*std::pow(-mypHE4HE4_PLCg_PIP2 + ypHE4HE4_PLCg_PIP2, 2)/std::pow(sigma_ypHE4HE4_PLCg_PIP2, 2);
            break;
        case 574:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHGF_Met_Met_PLCg_PIP2, 2)) + 0.5*std::pow(-mypHGF_Met_Met_PLCg_PIP2 + ypHGF_Met_Met_PLCg_PIP2, 2)/std::pow(sigma_ypHGF_Met_Met_PLCg_PIP2, 2);
            break;
        case 575:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHGF_Met_HGF_Met_PLCg_PIP2, 2)) + 0.5*std::pow(-mypHGF_Met_HGF_Met_PLCg_PIP2 + ypHGF_Met_HGF_Met_PLCg_PIP2, 2)/std::pow(sigma_ypHGF_Met_HGF_Met_PLCg_PIP2, 2);
            break;
        case 576:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypPPrPPr_PLCg_PIP2, 2)) + 0.5*std::pow(-mypPPrPPr_PLCg_PIP2 + ypPPrPPr_PLCg_PIP2, 2)/std::pow(sigma_ypPPrPPr_PLCg_PIP2, 2);
            break;
        case 577:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypPPrPr_PLCg_PIP2, 2)) + 0.5*std::pow(-mypPPrPr_PLCg_PIP2 + ypPPrPr_PLCg_PIP2, 2)/std::pow(sigma_ypPPrPr_PLCg_PIP2, 2);
            break;
        case 578:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypFFrFFr_PLCg_PIP2, 2)) + 0.5*std::pow(-mypFFrFFr_PLCg_PIP2 + ypFFrFFr_PLCg_PIP2, 2)/std::pow(sigma_ypFFrFFr_PLCg_PIP2, 2);
            break;
        case 579:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypFFrFr_PLCg_PIP2, 2)) + 0.5*std::pow(-mypFFrFr_PLCg_PIP2 + ypFFrFr_PLCg_PIP2, 2)/std::pow(sigma_ypFFrFr_PLCg_PIP2, 2);
            break;
        case 580:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypIIrIr_IRS_PLCg_PIP2, 2)) + 0.5*std::pow(-mypIIrIr_IRS_PLCg_PIP2 + ypIIrIr_IRS_PLCg_PIP2, 2)/std::pow(sigma_ypIIrIr_IRS_PLCg_PIP2, 2);
            break;
        case 581:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypINS_Isr_Isr_IRS_PLCg_PIP2, 2)) + 0.5*std::pow(-mypINS_Isr_Isr_IRS_PLCg_PIP2 + ypINS_Isr_Isr_IRS_PLCg_PIP2, 2)/std::pow(sigma_ypINS_Isr_Isr_IRS_PLCg_PIP2, 2);
            break;
        case 582:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypIIrIrI_IRS_PLCg_PIP2, 2)) + 0.5*std::pow(-mypIIrIrI_IRS_PLCg_PIP2 + ypIIrIrI_IRS_PLCg_PIP2, 2)/std::pow(sigma_ypIIrIrI_IRS_PLCg_PIP2, 2);
            break;
        case 583:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypINS_Isr_Isr_INS_IRS_PLCg_PIP2, 2)) + 0.5*std::pow(-mypINS_Isr_Isr_INS_IRS_PLCg_PIP2 + ypINS_Isr_Isr_INS_IRS_PLCg_PIP2, 2)/std::pow(sigma_ypINS_Isr_Isr_INS_IRS_PLCg_PIP2, 2);
            break;
        case 584:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E2_PI3K1_PIP2, 2)) + 0.5*std::pow(-mypEE1E2_PI3K1_PIP2 + ypEE1E2_PI3K1_PIP2, 2)/std::pow(sigma_ypEE1E2_PI3K1_PIP2, 2);
            break;
        case 585:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1Ev3_PI3K1_PIP2, 2)) + 0.5*std::pow(-mypEE1Ev3_PI3K1_PIP2 + ypEE1Ev3_PI3K1_PIP2, 2)/std::pow(sigma_ypEE1Ev3_PI3K1_PIP2, 2);
            break;
        case 586:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E1_PI3K1_PIP2, 2)) + 0.5*std::pow(-mypEE1E1_PI3K1_PIP2 + ypEE1E1_PI3K1_PIP2, 2)/std::pow(sigma_ypEE1E1_PI3K1_PIP2, 2);
            break;
        case 587:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1EE1_PI3K1_PIP2, 2)) + 0.5*std::pow(-mypEE1EE1_PI3K1_PIP2 + ypEE1EE1_PI3K1_PIP2, 2)/std::pow(sigma_ypEE1EE1_PI3K1_PIP2, 2);
            break;
        case 588:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E3_PI3K1_PIP2, 2)) + 0.5*std::pow(-mypEE1E3_PI3K1_PIP2 + ypEE1E3_PI3K1_PIP2, 2)/std::pow(sigma_ypEE1E3_PI3K1_PIP2, 2);
            break;
        case 589:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1HE3_PI3K1_PIP2, 2)) + 0.5*std::pow(-mypEE1HE3_PI3K1_PIP2 + ypEE1HE3_PI3K1_PIP2, 2)/std::pow(sigma_ypEE1HE3_PI3K1_PIP2, 2);
            break;
        case 590:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E4_PI3K1_PIP2, 2)) + 0.5*std::pow(-mypEE1E4_PI3K1_PIP2 + ypEE1E4_PI3K1_PIP2, 2)/std::pow(sigma_ypEE1E4_PI3K1_PIP2, 2);
            break;
        case 591:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1HE4_PI3K1_PIP2, 2)) + 0.5*std::pow(-mypEE1HE4_PI3K1_PIP2 + ypEE1HE4_PI3K1_PIP2, 2)/std::pow(sigma_ypEE1HE4_PI3K1_PIP2, 2);
            break;
        case 592:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE2HE3_PI3K1_PIP2, 2)) + 0.5*std::pow(-mypE2HE3_PI3K1_PIP2 + ypE2HE3_PI3K1_PIP2, 2)/std::pow(sigma_ypE2HE3_PI3K1_PIP2, 2);
            break;
        case 593:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE3Ev3_PI3K1_PIP2, 2)) + 0.5*std::pow(-mypHE3Ev3_PI3K1_PIP2 + ypHE3Ev3_PI3K1_PIP2, 2)/std::pow(sigma_ypHE3Ev3_PI3K1_PIP2, 2);
            break;
        case 594:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE1HE3_PI3K1_PIP2, 2)) + 0.5*std::pow(-mypE1HE3_PI3K1_PIP2 + ypE1HE3_PI3K1_PIP2, 2)/std::pow(sigma_ypE1HE3_PI3K1_PIP2, 2);
            break;
        case 595:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE3E4_PI3K1_PIP2, 2)) + 0.5*std::pow(-mypHE3E4_PI3K1_PIP2 + ypHE3E4_PI3K1_PIP2, 2)/std::pow(sigma_ypHE3E4_PI3K1_PIP2, 2);
            break;
        case 596:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE3HE4_PI3K1_PIP2, 2)) + 0.5*std::pow(-mypHE3HE4_PI3K1_PIP2 + ypHE3HE4_PI3K1_PIP2, 2)/std::pow(sigma_ypHE3HE4_PI3K1_PIP2, 2);
            break;
        case 597:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE2HE4_PI3K1_PIP2, 2)) + 0.5*std::pow(-mypE2HE4_PI3K1_PIP2 + ypE2HE4_PI3K1_PIP2, 2)/std::pow(sigma_ypE2HE4_PI3K1_PIP2, 2);
            break;
        case 598:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE4Ev3_PI3K1_PIP2, 2)) + 0.5*std::pow(-mypHE4Ev3_PI3K1_PIP2 + ypHE4Ev3_PI3K1_PIP2, 2)/std::pow(sigma_ypHE4Ev3_PI3K1_PIP2, 2);
            break;
        case 599:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE1HE4_PI3K1_PIP2, 2)) + 0.5*std::pow(-mypE1HE4_PI3K1_PIP2 + ypE1HE4_PI3K1_PIP2, 2)/std::pow(sigma_ypE1HE4_PI3K1_PIP2, 2);
            break;
        case 600:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE3HE4_PI3K1_PIP2, 2)) + 0.5*std::pow(-mypE3HE4_PI3K1_PIP2 + ypE3HE4_PI3K1_PIP2, 2)/std::pow(sigma_ypE3HE4_PI3K1_PIP2, 2);
            break;
        case 601:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE4E4_PI3K1_PIP2, 2)) + 0.5*std::pow(-mypHE4E4_PI3K1_PIP2 + ypHE4E4_PI3K1_PIP2, 2)/std::pow(sigma_ypHE4E4_PI3K1_PIP2, 2);
            break;
        case 602:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE4HE4_PI3K1_PIP2, 2)) + 0.5*std::pow(-mypHE4HE4_PI3K1_PIP2 + ypHE4HE4_PI3K1_PIP2, 2)/std::pow(sigma_ypHE4HE4_PI3K1_PIP2, 2);
            break;
        case 603:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHGF_Met_Met_PI3K1_PIP2, 2)) + 0.5*std::pow(-mypHGF_Met_Met_PI3K1_PIP2 + ypHGF_Met_Met_PI3K1_PIP2, 2)/std::pow(sigma_ypHGF_Met_Met_PI3K1_PIP2, 2);
            break;
        case 604:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHGF_Met_HGF_Met_PI3K1_PIP2, 2)) + 0.5*std::pow(-mypHGF_Met_HGF_Met_PI3K1_PIP2 + ypHGF_Met_HGF_Met_PI3K1_PIP2, 2)/std::pow(sigma_ypHGF_Met_HGF_Met_PI3K1_PIP2, 2);
            break;
        case 605:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypPPrPPr_PI3K1_PIP2, 2)) + 0.5*std::pow(-mypPPrPPr_PI3K1_PIP2 + ypPPrPPr_PI3K1_PIP2, 2)/std::pow(sigma_ypPPrPPr_PI3K1_PIP2, 2);
            break;
        case 606:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypPPrPr_PI3K1_PIP2, 2)) + 0.5*std::pow(-mypPPrPr_PI3K1_PIP2 + ypPPrPr_PI3K1_PIP2, 2)/std::pow(sigma_ypPPrPr_PI3K1_PIP2, 2);
            break;
        case 607:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypFFrFFr_PI3K1_PIP2, 2)) + 0.5*std::pow(-mypFFrFFr_PI3K1_PIP2 + ypFFrFFr_PI3K1_PIP2, 2)/std::pow(sigma_ypFFrFFr_PI3K1_PIP2, 2);
            break;
        case 608:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypFFrFr_PI3K1_PIP2, 2)) + 0.5*std::pow(-mypFFrFr_PI3K1_PIP2 + ypFFrFr_PI3K1_PIP2, 2)/std::pow(sigma_ypFFrFr_PI3K1_PIP2, 2);
            break;
        case 609:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypIIrIr_IRS_PI3K1_PIP2, 2)) + 0.5*std::pow(-mypIIrIr_IRS_PI3K1_PIP2 + ypIIrIr_IRS_PI3K1_PIP2, 2)/std::pow(sigma_ypIIrIr_IRS_PI3K1_PIP2, 2);
            break;
        case 610:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypINS_Isr_Isr_IRS_PI3K1_PIP2, 2)) + 0.5*std::pow(-mypINS_Isr_Isr_IRS_PI3K1_PIP2 + ypINS_Isr_Isr_IRS_PI3K1_PIP2, 2)/std::pow(sigma_ypINS_Isr_Isr_IRS_PI3K1_PIP2, 2);
            break;
        case 611:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypIIrIrI_IRS_PI3K1_PIP2, 2)) + 0.5*std::pow(-mypIIrIrI_IRS_PI3K1_PIP2 + ypIIrIrI_IRS_PI3K1_PIP2, 2)/std::pow(sigma_ypIIrIrI_IRS_PI3K1_PIP2, 2);
            break;
        case 612:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypINS_Isr_Isr_INS_IRS_PI3K1_PIP2, 2)) + 0.5*std::pow(-mypINS_Isr_Isr_INS_IRS_PI3K1_PIP2 + ypINS_Isr_Isr_INS_IRS_PI3K1_PIP2, 2)/std::pow(sigma_ypINS_Isr_Isr_INS_IRS_PI3K1_PIP2, 2);
            break;
        case 613:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E2_PI3K2_PIP, 2)) + 0.5*std::pow(-mypEE1E2_PI3K2_PIP + ypEE1E2_PI3K2_PIP, 2)/std::pow(sigma_ypEE1E2_PI3K2_PIP, 2);
            break;
        case 614:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1Ev3_PI3K2_PIP, 2)) + 0.5*std::pow(-mypEE1Ev3_PI3K2_PIP + ypEE1Ev3_PI3K2_PIP, 2)/std::pow(sigma_ypEE1Ev3_PI3K2_PIP, 2);
            break;
        case 615:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E1_PI3K2_PIP, 2)) + 0.5*std::pow(-mypEE1E1_PI3K2_PIP + ypEE1E1_PI3K2_PIP, 2)/std::pow(sigma_ypEE1E1_PI3K2_PIP, 2);
            break;
        case 616:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1EE1_PI3K2_PIP, 2)) + 0.5*std::pow(-mypEE1EE1_PI3K2_PIP + ypEE1EE1_PI3K2_PIP, 2)/std::pow(sigma_ypEE1EE1_PI3K2_PIP, 2);
            break;
        case 617:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E3_PI3K2_PIP, 2)) + 0.5*std::pow(-mypEE1E3_PI3K2_PIP + ypEE1E3_PI3K2_PIP, 2)/std::pow(sigma_ypEE1E3_PI3K2_PIP, 2);
            break;
        case 618:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1HE3_PI3K2_PIP, 2)) + 0.5*std::pow(-mypEE1HE3_PI3K2_PIP + ypEE1HE3_PI3K2_PIP, 2)/std::pow(sigma_ypEE1HE3_PI3K2_PIP, 2);
            break;
        case 619:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1E4_PI3K2_PIP, 2)) + 0.5*std::pow(-mypEE1E4_PI3K2_PIP + ypEE1E4_PI3K2_PIP, 2)/std::pow(sigma_ypEE1E4_PI3K2_PIP, 2);
            break;
        case 620:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEE1HE4_PI3K2_PIP, 2)) + 0.5*std::pow(-mypEE1HE4_PI3K2_PIP + ypEE1HE4_PI3K2_PIP, 2)/std::pow(sigma_ypEE1HE4_PI3K2_PIP, 2);
            break;
        case 621:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE2HE3_PI3K2_PIP, 2)) + 0.5*std::pow(-mypE2HE3_PI3K2_PIP + ypE2HE3_PI3K2_PIP, 2)/std::pow(sigma_ypE2HE3_PI3K2_PIP, 2);
            break;
        case 622:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE3Ev3_PI3K2_PIP, 2)) + 0.5*std::pow(-mypHE3Ev3_PI3K2_PIP + ypHE3Ev3_PI3K2_PIP, 2)/std::pow(sigma_ypHE3Ev3_PI3K2_PIP, 2);
            break;
        case 623:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE1HE3_PI3K2_PIP, 2)) + 0.5*std::pow(-mypE1HE3_PI3K2_PIP + ypE1HE3_PI3K2_PIP, 2)/std::pow(sigma_ypE1HE3_PI3K2_PIP, 2);
            break;
        case 624:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE3E4_PI3K2_PIP, 2)) + 0.5*std::pow(-mypHE3E4_PI3K2_PIP + ypHE3E4_PI3K2_PIP, 2)/std::pow(sigma_ypHE3E4_PI3K2_PIP, 2);
            break;
        case 625:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE3HE4_PI3K2_PIP, 2)) + 0.5*std::pow(-mypHE3HE4_PI3K2_PIP + ypHE3HE4_PI3K2_PIP, 2)/std::pow(sigma_ypHE3HE4_PI3K2_PIP, 2);
            break;
        case 626:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE2HE4_PI3K2_PIP, 2)) + 0.5*std::pow(-mypE2HE4_PI3K2_PIP + ypE2HE4_PI3K2_PIP, 2)/std::pow(sigma_ypE2HE4_PI3K2_PIP, 2);
            break;
        case 627:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE4Ev3_PI3K2_PIP, 2)) + 0.5*std::pow(-mypHE4Ev3_PI3K2_PIP + ypHE4Ev3_PI3K2_PIP, 2)/std::pow(sigma_ypHE4Ev3_PI3K2_PIP, 2);
            break;
        case 628:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE1HE4_PI3K2_PIP, 2)) + 0.5*std::pow(-mypE1HE4_PI3K2_PIP + ypE1HE4_PI3K2_PIP, 2)/std::pow(sigma_ypE1HE4_PI3K2_PIP, 2);
            break;
        case 629:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypE3HE4_PI3K2_PIP, 2)) + 0.5*std::pow(-mypE3HE4_PI3K2_PIP + ypE3HE4_PI3K2_PIP, 2)/std::pow(sigma_ypE3HE4_PI3K2_PIP, 2);
            break;
        case 630:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE4E4_PI3K2_PIP, 2)) + 0.5*std::pow(-mypHE4E4_PI3K2_PIP + ypHE4E4_PI3K2_PIP, 2)/std::pow(sigma_ypHE4E4_PI3K2_PIP, 2);
            break;
        case 631:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHE4HE4_PI3K2_PIP, 2)) + 0.5*std::pow(-mypHE4HE4_PI3K2_PIP + ypHE4HE4_PI3K2_PIP, 2)/std::pow(sigma_ypHE4HE4_PI3K2_PIP, 2);
            break;
        case 632:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHGF_Met_Met_PI3K2_PIP, 2)) + 0.5*std::pow(-mypHGF_Met_Met_PI3K2_PIP + ypHGF_Met_Met_PI3K2_PIP, 2)/std::pow(sigma_ypHGF_Met_Met_PI3K2_PIP, 2);
            break;
        case 633:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypHGF_Met_HGF_Met_PI3K2_PIP, 2)) + 0.5*std::pow(-mypHGF_Met_HGF_Met_PI3K2_PIP + ypHGF_Met_HGF_Met_PI3K2_PIP, 2)/std::pow(sigma_ypHGF_Met_HGF_Met_PI3K2_PIP, 2);
            break;
        case 634:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypPPrPPr_PI3K2_PIP, 2)) + 0.5*std::pow(-mypPPrPPr_PI3K2_PIP + ypPPrPPr_PI3K2_PIP, 2)/std::pow(sigma_ypPPrPPr_PI3K2_PIP, 2);
            break;
        case 635:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypPPrPr_PI3K2_PIP, 2)) + 0.5*std::pow(-mypPPrPr_PI3K2_PIP + ypPPrPr_PI3K2_PIP, 2)/std::pow(sigma_ypPPrPr_PI3K2_PIP, 2);
            break;
        case 636:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypFFrFFr_PI3K2_PIP, 2)) + 0.5*std::pow(-mypFFrFFr_PI3K2_PIP + ypFFrFFr_PI3K2_PIP, 2)/std::pow(sigma_ypFFrFFr_PI3K2_PIP, 2);
            break;
        case 637:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypFFrFr_PI3K2_PIP, 2)) + 0.5*std::pow(-mypFFrFr_PI3K2_PIP + ypFFrFr_PI3K2_PIP, 2)/std::pow(sigma_ypFFrFr_PI3K2_PIP, 2);
            break;
        case 638:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypIIrIr_IRS_PI3K2_PIP, 2)) + 0.5*std::pow(-mypIIrIr_IRS_PI3K2_PIP + ypIIrIr_IRS_PI3K2_PIP, 2)/std::pow(sigma_ypIIrIr_IRS_PI3K2_PIP, 2);
            break;
        case 639:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypINS_Isr_Isr_IRS_PI3K2_PIP, 2)) + 0.5*std::pow(-mypINS_Isr_Isr_IRS_PI3K2_PIP + ypINS_Isr_Isr_IRS_PI3K2_PIP, 2)/std::pow(sigma_ypINS_Isr_Isr_IRS_PI3K2_PIP, 2);
            break;
        case 640:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypIIrIrI_IRS_PI3K2_PIP, 2)) + 0.5*std::pow(-mypIIrIrI_IRS_PI3K2_PIP + ypIIrIrI_IRS_PI3K2_PIP, 2)/std::pow(sigma_ypIIrIrI_IRS_PI3K2_PIP, 2);
            break;
        case 641:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypINS_Isr_Isr_INS_IRS_PI3K2_PIP, 2)) + 0.5*std::pow(-mypINS_Isr_Isr_INS_IRS_PI3K2_PIP + ypINS_Isr_Isr_INS_IRS_PI3K2_PIP, 2)/std::pow(sigma_ypINS_Isr_Isr_INS_IRS_PI3K2_PIP, 2);
            break;
        case 642:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yIRS, 2)) + 0.5*std::pow(-myIRS + yIRS, 2)/std::pow(sigma_yIRS, 2);
            break;
        case 643:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ySp, 2)) + 0.5*std::pow(-mySp + ySp, 2)/std::pow(sigma_ySp, 2);
            break;
        case 644:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yCbl, 2)) + 0.5*std::pow(-myCbl + yCbl, 2)/std::pow(sigma_yCbl, 2);
            break;
        case 645:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yG2, 2)) + 0.5*std::pow(-myG2 + yG2, 2)/std::pow(sigma_yG2, 2);
            break;
        case 646:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yG2_SOS, 2)) + 0.5*std::pow(-myG2_SOS + yG2_SOS, 2)/std::pow(sigma_yG2_SOS, 2);
            break;
        case 647:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yG2_pSOS, 2)) + 0.5*std::pow(-myG2_pSOS + yG2_pSOS, 2)/std::pow(sigma_yG2_pSOS, 2);
            break;
        case 648:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yPLCg, 2)) + 0.5*std::pow(-myPLCg + yPLCg, 2)/std::pow(sigma_yPLCg, 2);
            break;
        case 649:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yPI3KC1, 2)) + 0.5*std::pow(-myPI3KC1 + yPI3KC1, 2)/std::pow(sigma_yPI3KC1, 2);
            break;
        case 650:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yPI3KR1, 2)) + 0.5*std::pow(-myPI3KR1 + yPI3KR1, 2)/std::pow(sigma_yPI3KR1, 2);
            break;
        case 651:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yPI3K1, 2)) + 0.5*std::pow(-myPI3K1 + yPI3K1, 2)/std::pow(sigma_yPI3K1, 2);
            break;
        case 652:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypPI3K1, 2)) + 0.5*std::pow(-mypPI3K1 + ypPI3K1, 2)/std::pow(sigma_ypPI3K1, 2);
            break;
        case 653:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yPI3K2, 2)) + 0.5*std::pow(-myPI3K2 + yPI3K2, 2)/std::pow(sigma_yPI3K2, 2);
            break;
        case 654:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ymTORC1, 2)) + 0.5*std::pow(-mymTORC1 + ymTORC1, 2)/std::pow(sigma_ymTORC1, 2);
            break;
        case 655:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ymTORC1active, 2)) + 0.5*std::pow(-mymTORC1active + ymTORC1active, 2)/std::pow(sigma_ymTORC1active, 2);
            break;
        case 656:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yPIP, 2)) + 0.5*std::pow(-myPIP + yPIP, 2)/std::pow(sigma_yPIP, 2);
            break;
        case 657:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yPI3P, 2)) + 0.5*std::pow(-myPI3P + yPI3P, 2)/std::pow(sigma_yPI3P, 2);
            break;
        case 658:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yDAG, 2)) + 0.5*std::pow(-myDAG + yDAG, 2)/std::pow(sigma_yDAG, 2);
            break;
        case 659:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yGRP, 2)) + 0.5*std::pow(-myGRP + yGRP, 2)/std::pow(sigma_yGRP, 2);
            break;
        case 660:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yDAG_GRP, 2)) + 0.5*std::pow(-myDAG_GRP + yDAG_GRP, 2)/std::pow(sigma_yDAG_GRP, 2);
            break;
        case 661:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yRasT, 2)) + 0.5*std::pow(-myRasT + yRasT, 2)/std::pow(sigma_yRasT, 2);
            break;
        case 662:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yRasD, 2)) + 0.5*std::pow(-myRasD + yRasD, 2)/std::pow(sigma_yRasD, 2);
            break;
        case 663:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yNF1, 2)) + 0.5*std::pow(-myNF1 + yNF1, 2)/std::pow(sigma_yNF1, 2);
            break;
        case 664:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypNF1, 2)) + 0.5*std::pow(-mypNF1 + ypNF1, 2)/std::pow(sigma_ypNF1, 2);
            break;
        case 665:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypCRaf, 2)) + 0.5*std::pow(-mypCRaf + ypCRaf, 2)/std::pow(sigma_ypCRaf, 2);
            break;
        case 666:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yCRaf, 2)) + 0.5*std::pow(-myCRaf + yCRaf, 2)/std::pow(sigma_yCRaf, 2);
            break;
        case 667:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yRasT_CRaf, 2)) + 0.5*std::pow(-myRasT_CRaf + yRasT_CRaf, 2)/std::pow(sigma_yRasT_CRaf, 2);
            break;
        case 668:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yBRaf, 2)) + 0.5*std::pow(-myBRaf + yBRaf, 2)/std::pow(sigma_yBRaf, 2);
            break;
        case 669:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yRasT_CRaf_BRaf, 2)) + 0.5*std::pow(-myRasT_CRaf_BRaf + yRasT_CRaf_BRaf, 2)/std::pow(sigma_yRasT_CRaf_BRaf, 2);
            break;
        case 670:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMEK, 2)) + 0.5*std::pow(-myMEK + yMEK, 2)/std::pow(sigma_yMEK, 2);
            break;
        case 671:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypMEK, 2)) + 0.5*std::pow(-mypMEK + ypMEK, 2)/std::pow(sigma_ypMEK, 2);
            break;
        case 672:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yppMEK, 2)) + 0.5*std::pow(-myppMEK + yppMEK, 2)/std::pow(sigma_yppMEK, 2);
            break;
        case 673:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMKP3, 2)) + 0.5*std::pow(-myMKP3 + yMKP3, 2)/std::pow(sigma_yMKP3, 2);
            break;
        case 674:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yERKnuc, 2)) + 0.5*std::pow(-myERKnuc + yERKnuc, 2)/std::pow(sigma_yERKnuc, 2);
            break;
        case 675:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yppERKnuc, 2)) + 0.5*std::pow(-myppERKnuc + yppERKnuc, 2)/std::pow(sigma_yppERKnuc, 2);
            break;
        case 676:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yRSK, 2)) + 0.5*std::pow(-myRSK + yRSK, 2)/std::pow(sigma_yRSK, 2);
            break;
        case 677:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypRSK, 2)) + 0.5*std::pow(-mypRSK + ypRSK, 2)/std::pow(sigma_ypRSK, 2);
            break;
        case 678:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypRSKnuc, 2)) + 0.5*std::pow(-mypRSKnuc + ypRSKnuc, 2)/std::pow(sigma_ypRSKnuc, 2);
            break;
        case 679:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMKP1, 2)) + 0.5*std::pow(-myMKP1 + yMKP1, 2)/std::pow(sigma_yMKP1, 2);
            break;
        case 680:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypMKP1, 2)) + 0.5*std::pow(-mypMKP1 + ypMKP1, 2)/std::pow(sigma_ypMKP1, 2);
            break;
        case 681:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ycFos, 2)) + 0.5*std::pow(-mycFos + ycFos, 2)/std::pow(sigma_ycFos, 2);
            break;
        case 682:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypcFos, 2)) + 0.5*std::pow(-mypcFos + ypcFos, 2)/std::pow(sigma_ypcFos, 2);
            break;
        case 683:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ycJun, 2)) + 0.5*std::pow(-mycJun + ycJun, 2)/std::pow(sigma_ycJun, 2);
            break;
        case 684:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypcFos_cJun, 2)) + 0.5*std::pow(-mypcFos_cJun + ypcFos_cJun, 2)/std::pow(sigma_ypcFos_cJun, 2);
            break;
        case 685:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ycMyc, 2)) + 0.5*std::pow(-mycMyc + ycMyc, 2)/std::pow(sigma_ycMyc, 2);
            break;
        case 686:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ybCATENINnuc, 2)) + 0.5*std::pow(-mybCATENINnuc + ybCATENINnuc, 2)/std::pow(sigma_ybCATENINnuc, 2);
            break;
        case 687:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ybCATENIN, 2)) + 0.5*std::pow(-mybCATENIN + ybCATENIN, 2)/std::pow(sigma_ybCATENIN, 2);
            break;
        case 688:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypbCATENIN, 2)) + 0.5*std::pow(-mypbCATENIN + ypbCATENIN, 2)/std::pow(sigma_ypbCATENIN, 2);
            break;
        case 689:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yIP3, 2)) + 0.5*std::pow(-myIP3 + yIP3, 2)/std::pow(sigma_yIP3, 2);
            break;
        case 690:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yPIP2, 2)) + 0.5*std::pow(-myPIP2 + yPIP2, 2)/std::pow(sigma_yPIP2, 2);
            break;
        case 691:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yPIP3, 2)) + 0.5*std::pow(-myPIP3 + yPIP3, 2)/std::pow(sigma_yPIP3, 2);
            break;
        case 692:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yPTEN, 2)) + 0.5*std::pow(-myPTEN + yPTEN, 2)/std::pow(sigma_yPTEN, 2);
            break;
        case 693:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yPIP3_AKT, 2)) + 0.5*std::pow(-myPIP3_AKT + yPIP3_AKT, 2)/std::pow(sigma_yPIP3_AKT, 2);
            break;
        case 694:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yAKT, 2)) + 0.5*std::pow(-myAKT + yAKT, 2)/std::pow(sigma_yAKT, 2);
            break;
        case 695:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypAKT, 2)) + 0.5*std::pow(-mypAKT + ypAKT, 2)/std::pow(sigma_ypAKT, 2);
            break;
        case 696:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yppAKT, 2)) + 0.5*std::pow(-myppAKT + yppAKT, 2)/std::pow(sigma_yppAKT, 2);
            break;
        case 697:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yPDK1, 2)) + 0.5*std::pow(-myPDK1 + yPDK1, 2)/std::pow(sigma_yPDK1, 2);
            break;
        case 698:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yPIP3_PDK1, 2)) + 0.5*std::pow(-myPIP3_PDK1 + yPIP3_PDK1, 2)/std::pow(sigma_yPIP3_PDK1, 2);
            break;
        case 699:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yPIP3_pAKT, 2)) + 0.5*std::pow(-myPIP3_pAKT + yPIP3_pAKT, 2)/std::pow(sigma_yPIP3_pAKT, 2);
            break;
        case 700:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yRictor, 2)) + 0.5*std::pow(-myRictor + yRictor, 2)/std::pow(sigma_yRictor, 2);
            break;
        case 701:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ymTOR, 2)) + 0.5*std::pow(-mymTOR + ymTOR, 2)/std::pow(sigma_ymTOR, 2);
            break;
        case 702:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ymTORC2, 2)) + 0.5*std::pow(-mymTORC2 + ymTORC2, 2)/std::pow(sigma_ymTORC2, 2);
            break;
        case 703:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yPIP3_ppAKT, 2)) + 0.5*std::pow(-myPIP3_ppAKT + yPIP3_ppAKT, 2)/std::pow(sigma_yPIP3_ppAKT, 2);
            break;
        case 704:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yGSK3b, 2)) + 0.5*std::pow(-myGSK3b + yGSK3b, 2)/std::pow(sigma_yGSK3b, 2);
            break;
        case 705:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypGSK3b, 2)) + 0.5*std::pow(-mypGSK3b + ypGSK3b, 2)/std::pow(sigma_ypGSK3b, 2);
            break;
        case 706:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yTSC1, 2)) + 0.5*std::pow(-myTSC1 + yTSC1, 2)/std::pow(sigma_yTSC1, 2);
            break;
        case 707:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yTSC2, 2)) + 0.5*std::pow(-myTSC2 + yTSC2, 2)/std::pow(sigma_yTSC2, 2);
            break;
        case 708:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypTSC2, 2)) + 0.5*std::pow(-mypTSC2 + ypTSC2, 2)/std::pow(sigma_ypTSC2, 2);
            break;
        case 709:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yTSC, 2)) + 0.5*std::pow(-myTSC + yTSC, 2)/std::pow(sigma_yTSC, 2);
            break;
        case 710:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yPKC, 2)) + 0.5*std::pow(-myPKC + yPKC, 2)/std::pow(sigma_yPKC, 2);
            break;
        case 711:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yDAG_PKC, 2)) + 0.5*std::pow(-myDAG_PKC + yDAG_PKC, 2)/std::pow(sigma_yDAG_PKC, 2);
            break;
        case 712:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypRKIP, 2)) + 0.5*std::pow(-mypRKIP + ypRKIP, 2)/std::pow(sigma_ypRKIP, 2);
            break;
        case 713:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yRKIP, 2)) + 0.5*std::pow(-myRKIP + yRKIP, 2)/std::pow(sigma_yRKIP, 2);
            break;
        case 714:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yRKIP_CRaf, 2)) + 0.5*std::pow(-myRKIP_CRaf + yRKIP_CRaf, 2)/std::pow(sigma_yRKIP_CRaf, 2);
            break;
        case 715:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yERK, 2)) + 0.5*std::pow(-myERK + yERK, 2)/std::pow(sigma_yERK, 2);
            break;
        case 716:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypERK, 2)) + 0.5*std::pow(-mypERK + ypERK, 2)/std::pow(sigma_ypERK, 2);
            break;
        case 717:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yppERK, 2)) + 0.5*std::pow(-myppERK + yppERK, 2)/std::pow(sigma_yppERK, 2);
            break;
        case 718:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yFOXO, 2)) + 0.5*std::pow(-myFOXO + yFOXO, 2)/std::pow(sigma_yFOXO, 2);
            break;
        case 719:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypFOXO, 2)) + 0.5*std::pow(-mypFOXO + ypFOXO, 2)/std::pow(sigma_ypFOXO, 2);
            break;
        case 720:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yRhebD, 2)) + 0.5*std::pow(-myRhebD + yRhebD, 2)/std::pow(sigma_yRhebD, 2);
            break;
        case 721:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yRhebT, 2)) + 0.5*std::pow(-myRhebT + yRhebT, 2)/std::pow(sigma_yRhebT, 2);
            break;
        case 722:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yRaptor, 2)) + 0.5*std::pow(-myRaptor + yRaptor, 2)/std::pow(sigma_yRaptor, 2);
            break;
        case 723:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yS6K, 2)) + 0.5*std::pow(-myS6K + yS6K, 2)/std::pow(sigma_yS6K, 2);
            break;
        case 724:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypS6K, 2)) + 0.5*std::pow(-mypS6K + ypS6K, 2)/std::pow(sigma_ypS6K, 2);
            break;
        case 725:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yEIF4EBP1, 2)) + 0.5*std::pow(-myEIF4EBP1 + yEIF4EBP1, 2)/std::pow(sigma_yEIF4EBP1, 2);
            break;
        case 726:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypEIF4EBP1, 2)) + 0.5*std::pow(-mypEIF4EBP1 + ypEIF4EBP1, 2)/std::pow(sigma_ypEIF4EBP1, 2);
            break;
        case 727:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ySOS, 2)) + 0.5*std::pow(-mySOS + ySOS, 2)/std::pow(sigma_ySOS, 2);
            break;
        case 728:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yG2_SOS_ppERK, 2)) + 0.5*std::pow(-myG2_SOS_ppERK + yG2_SOS_ppERK, 2)/std::pow(sigma_yG2_SOS_ppERK, 2);
            break;
        case 729:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yCRaf_ppERK, 2)) + 0.5*std::pow(-myCRaf_ppERK + yCRaf_ppERK, 2)/std::pow(sigma_yCRaf_ppERK, 2);
            break;
        case 730:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yRasD_DAG_GRP, 2)) + 0.5*std::pow(-myRasD_DAG_GRP + yRasD_DAG_GRP, 2)/std::pow(sigma_yRasD_DAG_GRP, 2);
            break;
        case 731:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yRasT_NF1, 2)) + 0.5*std::pow(-myRasT_NF1 + yRasT_NF1, 2)/std::pow(sigma_yRasT_NF1, 2);
            break;
        case 732:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yNF1_ppERK, 2)) + 0.5*std::pow(-myNF1_ppERK + yNF1_ppERK, 2)/std::pow(sigma_yNF1_ppERK, 2);
            break;
        case 733:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMEK_RasT_CRaf_BRaf, 2)) + 0.5*std::pow(-myMEK_RasT_CRaf_BRaf + yMEK_RasT_CRaf_BRaf, 2)/std::pow(sigma_yMEK_RasT_CRaf_BRaf, 2);
            break;
        case 734:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypMEK_RasT_CRaf_BRaf, 2)) + 0.5*std::pow(-mypMEK_RasT_CRaf_BRaf + ypMEK_RasT_CRaf_BRaf, 2)/std::pow(sigma_ypMEK_RasT_CRaf_BRaf, 2);
            break;
        case 735:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yERK_ppMEK, 2)) + 0.5*std::pow(-myERK_ppMEK + yERK_ppMEK, 2)/std::pow(sigma_yERK_ppMEK, 2);
            break;
        case 736:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypERK_ppMEK, 2)) + 0.5*std::pow(-mypERK_ppMEK + ypERK_ppMEK, 2)/std::pow(sigma_ypERK_ppMEK, 2);
            break;
        case 737:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yRSK_ppERK, 2)) + 0.5*std::pow(-myRSK_ppERK + yRSK_ppERK, 2)/std::pow(sigma_yRSK_ppERK, 2);
            break;
        case 738:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypRSKnuc_MKP1, 2)) + 0.5*std::pow(-mypRSKnuc_MKP1 + ypRSKnuc_MKP1, 2)/std::pow(sigma_ypRSKnuc_MKP1, 2);
            break;
        case 739:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yppERKnuc_MKP1, 2)) + 0.5*std::pow(-myppERKnuc_MKP1 + yppERKnuc_MKP1, 2)/std::pow(sigma_yppERKnuc_MKP1, 2);
            break;
        case 740:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ycFos_pRSKnuc, 2)) + 0.5*std::pow(-mycFos_pRSKnuc + ycFos_pRSKnuc, 2)/std::pow(sigma_ycFos_pRSKnuc, 2);
            break;
        case 741:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ycFos_ppERKnuc, 2)) + 0.5*std::pow(-mycFos_ppERKnuc + ycFos_ppERKnuc, 2)/std::pow(sigma_ycFos_ppERKnuc, 2);
            break;
        case 742:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yRKIP_DAG_PKC, 2)) + 0.5*std::pow(-myRKIP_DAG_PKC + yRKIP_DAG_PKC, 2)/std::pow(sigma_yRKIP_DAG_PKC, 2);
            break;
        case 743:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yPIP3_PTEN, 2)) + 0.5*std::pow(-myPIP3_PTEN + yPIP3_PTEN, 2)/std::pow(sigma_yPIP3_PTEN, 2);
            break;
        case 744:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yPIP3_AKT_PIP3_PDK1, 2)) + 0.5*std::pow(-myPIP3_AKT_PIP3_PDK1 + yPIP3_AKT_PIP3_PDK1, 2)/std::pow(sigma_yPIP3_AKT_PIP3_PDK1, 2);
            break;
        case 745:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yPIP3_pAKT_mTORC2, 2)) + 0.5*std::pow(-myPIP3_pAKT_mTORC2 + yPIP3_pAKT_mTORC2, 2)/std::pow(sigma_yPIP3_pAKT_mTORC2, 2);
            break;
        case 746:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yGSK3b_ppAKT, 2)) + 0.5*std::pow(-myGSK3b_ppAKT + yGSK3b_ppAKT, 2)/std::pow(sigma_yGSK3b_ppAKT, 2);
            break;
        case 747:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yTSC2_ppAKT, 2)) + 0.5*std::pow(-myTSC2_ppAKT + yTSC2_ppAKT, 2)/std::pow(sigma_yTSC2_ppAKT, 2);
            break;
        case 748:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yTSC2_ppERK, 2)) + 0.5*std::pow(-myTSC2_ppERK + yTSC2_ppERK, 2)/std::pow(sigma_yTSC2_ppERK, 2);
            break;
        case 749:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yRhebT_TSC, 2)) + 0.5*std::pow(-myRhebT_TSC + yRhebT_TSC, 2)/std::pow(sigma_yRhebT_TSC, 2);
            break;
        case 750:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yEIF4EBP1_mTORC1active, 2)) + 0.5*std::pow(-myEIF4EBP1_mTORC1active + yEIF4EBP1_mTORC1active, 2)/std::pow(sigma_yEIF4EBP1_mTORC1active, 2);
            break;
        case 751:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yS6K_mTORC1active, 2)) + 0.5*std::pow(-myS6K_mTORC1active + yS6K_mTORC1active, 2)/std::pow(sigma_yS6K_mTORC1active, 2);
            break;
        case 752:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yFOXO_ppAKT, 2)) + 0.5*std::pow(-myFOXO_ppAKT + yFOXO_ppAKT, 2)/std::pow(sigma_yFOXO_ppAKT, 2);
            break;
        case 753:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yPI3K1_mTORC1active, 2)) + 0.5*std::pow(-myPI3K1_mTORC1active + yPI3K1_mTORC1active, 2)/std::pow(sigma_yPI3K1_mTORC1active, 2);
            break;
        case 754:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypERK_MKP3, 2)) + 0.5*std::pow(-mypERK_MKP3 + ypERK_MKP3, 2)/std::pow(sigma_ypERK_MKP3, 2);
            break;
        case 755:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yppERK_MKP3, 2)) + 0.5*std::pow(-myppERK_MKP3 + yppERK_MKP3, 2)/std::pow(sigma_yppERK_MKP3, 2);
            break;
        case 756:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yppERKnuc_pMKP1, 2)) + 0.5*std::pow(-myppERKnuc_pMKP1 + yppERKnuc_pMKP1, 2)/std::pow(sigma_yppERKnuc_pMKP1, 2);
            break;
        case 757:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yRasT_BRaf, 2)) + 0.5*std::pow(-myRasT_BRaf + yRasT_BRaf, 2)/std::pow(sigma_yRasT_BRaf, 2);
            break;
        case 758:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yRasT_BRaf_BRaf, 2)) + 0.5*std::pow(-myRasT_BRaf_BRaf + yRasT_BRaf_BRaf, 2)/std::pow(sigma_yRasT_BRaf_BRaf, 2);
            break;
        case 759:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMEK_RasT_BRaf_BRaf, 2)) + 0.5*std::pow(-myMEK_RasT_BRaf_BRaf + yMEK_RasT_BRaf_BRaf, 2)/std::pow(sigma_yMEK_RasT_BRaf_BRaf, 2);
            break;
        case 760:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypMEK_RasT_BRaf_BRaf, 2)) + 0.5*std::pow(-mypMEK_RasT_BRaf_BRaf + ypMEK_RasT_BRaf_BRaf, 2)/std::pow(sigma_ypMEK_RasT_BRaf_BRaf, 2);
            break;
        case 761:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yEIF4E, 2)) + 0.5*std::pow(-myEIF4E + yEIF4E, 2)/std::pow(sigma_yEIF4E, 2);
            break;
        case 762:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yEIF4EBP1_EIF4E, 2)) + 0.5*std::pow(-myEIF4EBP1_EIF4E + yEIF4EBP1_EIF4E, 2)/std::pow(sigma_yEIF4EBP1_EIF4E, 2);
            break;
        case 763:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yRasT_CRaf_CRaf, 2)) + 0.5*std::pow(-myRasT_CRaf_CRaf + yRasT_CRaf_CRaf, 2)/std::pow(sigma_yRasT_CRaf_CRaf, 2);
            break;
        case 764:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMEK_RasT_CRaf_CRaf, 2)) + 0.5*std::pow(-myMEK_RasT_CRaf_CRaf + yMEK_RasT_CRaf_CRaf, 2)/std::pow(sigma_yMEK_RasT_CRaf_CRaf, 2);
            break;
        case 765:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ypMEK_RasT_CRaf_CRaf, 2)) + 0.5*std::pow(-mypMEK_RasT_CRaf_CRaf + ypMEK_RasT_CRaf_CRaf, 2)/std::pow(sigma_ypMEK_RasT_CRaf_CRaf, 2);
            break;
        case 766:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yFOXOnuc, 2)) + 0.5*std::pow(-myFOXOnuc + yFOXOnuc, 2)/std::pow(sigma_yFOXOnuc, 2);
            break;
        case 767:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMEKi, 2)) + 0.5*std::pow(-myMEKi + yMEKi, 2)/std::pow(sigma_yMEKi, 2);
            break;
        case 768:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMEKi_ppMEK, 2)) + 0.5*std::pow(-myMEKi_ppMEK + yMEKi_ppMEK, 2)/std::pow(sigma_yMEKi_ppMEK, 2);
            break;
        case 769:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yAKTi, 2)) + 0.5*std::pow(-myAKTi + yAKTi, 2)/std::pow(sigma_yAKTi, 2);
            break;
        case 770:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yAKTi_AKT, 2)) + 0.5*std::pow(-myAKTi_AKT + yAKTi_AKT, 2)/std::pow(sigma_yAKTi_AKT, 2);
            break;
        case 771:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ymT, 2)) + 0.5*std::pow(-mymT + ymT, 2)/std::pow(sigma_ymT, 2);
            break;
        case 772:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yEIF4E_mT, 2)) + 0.5*std::pow(-myEIF4E_mT + yEIF4E_mT, 2)/std::pow(sigma_yEIF4E_mT, 2);
            break;
        case 773:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_TP53, 2)) + 0.5*std::pow(-mym_TP53 + ym_TP53, 2)/std::pow(sigma_ym_TP53, 2);
            break;
        case 774:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_MDM2, 2)) + 0.5*std::pow(-mym_MDM2 + ym_MDM2, 2)/std::pow(sigma_ym_MDM2, 2);
            break;
        case 775:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_PPM1D, 2)) + 0.5*std::pow(-mym_PPM1D + ym_PPM1D, 2)/std::pow(sigma_ym_PPM1D, 2);
            break;
        case 776:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_ATM, 2)) + 0.5*std::pow(-mym_ATM + ym_ATM, 2)/std::pow(sigma_ym_ATM, 2);
            break;
        case 777:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_ATR, 2)) + 0.5*std::pow(-mym_ATR + ym_ATR, 2)/std::pow(sigma_ym_ATR, 2);
            break;
        case 778:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_RB1, 2)) + 0.5*std::pow(-mym_RB1 + ym_RB1, 2)/std::pow(sigma_ym_RB1, 2);
            break;
        case 779:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_E2F1, 2)) + 0.5*std::pow(-mym_E2F1 + ym_E2F1, 2)/std::pow(sigma_ym_E2F1, 2);
            break;
        case 780:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_E2F2, 2)) + 0.5*std::pow(-mym_E2F2 + ym_E2F2, 2)/std::pow(sigma_ym_E2F2, 2);
            break;
        case 781:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_E2F3, 2)) + 0.5*std::pow(-mym_E2F3 + ym_E2F3, 2)/std::pow(sigma_ym_E2F3, 2);
            break;
        case 782:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_CCND1, 2)) + 0.5*std::pow(-mym_CCND1 + ym_CCND1, 2)/std::pow(sigma_ym_CCND1, 2);
            break;
        case 783:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_CCND2, 2)) + 0.5*std::pow(-mym_CCND2 + ym_CCND2, 2)/std::pow(sigma_ym_CCND2, 2);
            break;
        case 784:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_CCND3, 2)) + 0.5*std::pow(-mym_CCND3 + ym_CCND3, 2)/std::pow(sigma_ym_CCND3, 2);
            break;
        case 785:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_CCNE1, 2)) + 0.5*std::pow(-mym_CCNE1 + ym_CCNE1, 2)/std::pow(sigma_ym_CCNE1, 2);
            break;
        case 786:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_CCNE2, 2)) + 0.5*std::pow(-mym_CCNE2 + ym_CCNE2, 2)/std::pow(sigma_ym_CCNE2, 2);
            break;
        case 787:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_SKP2, 2)) + 0.5*std::pow(-mym_SKP2 + ym_SKP2, 2)/std::pow(sigma_ym_SKP2, 2);
            break;
        case 788:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_CDC25A, 2)) + 0.5*std::pow(-mym_CDC25A + ym_CDC25A, 2)/std::pow(sigma_ym_CDC25A, 2);
            break;
        case 789:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_CDC25B, 2)) + 0.5*std::pow(-mym_CDC25B + ym_CDC25B, 2)/std::pow(sigma_ym_CDC25B, 2);
            break;
        case 790:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_CDC25C, 2)) + 0.5*std::pow(-mym_CDC25C + ym_CDC25C, 2)/std::pow(sigma_ym_CDC25C, 2);
            break;
        case 791:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_CCNA2, 2)) + 0.5*std::pow(-mym_CCNA2 + ym_CCNA2, 2)/std::pow(sigma_ym_CCNA2, 2);
            break;
        case 792:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_CDKN1B, 2)) + 0.5*std::pow(-mym_CDKN1B + ym_CDKN1B, 2)/std::pow(sigma_ym_CDKN1B, 2);
            break;
        case 793:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_CDH1, 2)) + 0.5*std::pow(-mym_CDH1 + ym_CDH1, 2)/std::pow(sigma_ym_CDH1, 2);
            break;
        case 794:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_CCNB1, 2)) + 0.5*std::pow(-mym_CCNB1 + ym_CCNB1, 2)/std::pow(sigma_ym_CCNB1, 2);
            break;
        case 795:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_CDC20, 2)) + 0.5*std::pow(-mym_CDC20 + ym_CDC20, 2)/std::pow(sigma_ym_CDC20, 2);
            break;
        case 796:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_WEE1, 2)) + 0.5*std::pow(-mym_WEE1 + ym_WEE1, 2)/std::pow(sigma_ym_WEE1, 2);
            break;
        case 797:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_CHEK1, 2)) + 0.5*std::pow(-mym_CHEK1 + ym_CHEK1, 2)/std::pow(sigma_ym_CHEK1, 2);
            break;
        case 798:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_CDKN1A, 2)) + 0.5*std::pow(-mym_CDKN1A + ym_CDKN1A, 2)/std::pow(sigma_ym_CDKN1A, 2);
            break;
        case 799:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_CDK1, 2)) + 0.5*std::pow(-mym_CDK1 + ym_CDK1, 2)/std::pow(sigma_ym_CDK1, 2);
            break;
        case 800:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_CDK2, 2)) + 0.5*std::pow(-mym_CDK2 + ym_CDK2, 2)/std::pow(sigma_ym_CDK2, 2);
            break;
        case 801:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_CDK4, 2)) + 0.5*std::pow(-mym_CDK4 + ym_CDK4, 2)/std::pow(sigma_ym_CDK4, 2);
            break;
        case 802:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_CDK6, 2)) + 0.5*std::pow(-mym_CDK6 + ym_CDK6, 2)/std::pow(sigma_ym_CDK6, 2);
            break;
        case 803:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_TNFSF10, 2)) + 0.5*std::pow(-mym_TNFSF10 + ym_TNFSF10, 2)/std::pow(sigma_ym_TNFSF10, 2);
            break;
        case 804:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_TNFRSF10A, 2)) + 0.5*std::pow(-mym_TNFRSF10A + ym_TNFRSF10A, 2)/std::pow(sigma_ym_TNFRSF10A, 2);
            break;
        case 805:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_TNFRSF10B, 2)) + 0.5*std::pow(-mym_TNFRSF10B + ym_TNFRSF10B, 2)/std::pow(sigma_ym_TNFRSF10B, 2);
            break;
        case 806:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_CFLAR, 2)) + 0.5*std::pow(-mym_CFLAR + ym_CFLAR, 2)/std::pow(sigma_ym_CFLAR, 2);
            break;
        case 807:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_CASP8, 2)) + 0.5*std::pow(-mym_CASP8 + ym_CASP8, 2)/std::pow(sigma_ym_CASP8, 2);
            break;
        case 808:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_CASP10, 2)) + 0.5*std::pow(-mym_CASP10 + ym_CASP10, 2)/std::pow(sigma_ym_CASP10, 2);
            break;
        case 809:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_BFAR, 2)) + 0.5*std::pow(-mym_BFAR + ym_BFAR, 2)/std::pow(sigma_ym_BFAR, 2);
            break;
        case 810:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_CASP3, 2)) + 0.5*std::pow(-mym_CASP3 + ym_CASP3, 2)/std::pow(sigma_ym_CASP3, 2);
            break;
        case 811:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_CASP7, 2)) + 0.5*std::pow(-mym_CASP7 + ym_CASP7, 2)/std::pow(sigma_ym_CASP7, 2);
            break;
        case 812:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_CASP6, 2)) + 0.5*std::pow(-mym_CASP6 + ym_CASP6, 2)/std::pow(sigma_ym_CASP6, 2);
            break;
        case 813:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_XIAP, 2)) + 0.5*std::pow(-mym_XIAP + ym_XIAP, 2)/std::pow(sigma_ym_XIAP, 2);
            break;
        case 814:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_PARP1, 2)) + 0.5*std::pow(-mym_PARP1 + ym_PARP1, 2)/std::pow(sigma_ym_PARP1, 2);
            break;
        case 815:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_BID, 2)) + 0.5*std::pow(-mym_BID + ym_BID, 2)/std::pow(sigma_ym_BID, 2);
            break;
        case 816:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_BCL2, 2)) + 0.5*std::pow(-mym_BCL2 + ym_BCL2, 2)/std::pow(sigma_ym_BCL2, 2);
            break;
        case 817:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_BCL2L1, 2)) + 0.5*std::pow(-mym_BCL2L1 + ym_BCL2L1, 2)/std::pow(sigma_ym_BCL2L1, 2);
            break;
        case 818:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_MCL1, 2)) + 0.5*std::pow(-mym_MCL1 + ym_MCL1, 2)/std::pow(sigma_ym_MCL1, 2);
            break;
        case 819:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_BAX, 2)) + 0.5*std::pow(-mym_BAX + ym_BAX, 2)/std::pow(sigma_ym_BAX, 2);
            break;
        case 820:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_CYCS, 2)) + 0.5*std::pow(-mym_CYCS + ym_CYCS, 2)/std::pow(sigma_ym_CYCS, 2);
            break;
        case 821:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_DIABLO, 2)) + 0.5*std::pow(-mym_DIABLO + ym_DIABLO, 2)/std::pow(sigma_ym_DIABLO, 2);
            break;
        case 822:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_APAF1, 2)) + 0.5*std::pow(-mym_APAF1 + ym_APAF1, 2)/std::pow(sigma_ym_APAF1, 2);
            break;
        case 823:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_CASP9, 2)) + 0.5*std::pow(-mym_CASP9 + ym_CASP9, 2)/std::pow(sigma_ym_CASP9, 2);
            break;
        case 824:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_BAD, 2)) + 0.5*std::pow(-mym_BAD + ym_BAD, 2)/std::pow(sigma_ym_BAD, 2);
            break;
        case 825:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_BBC3, 2)) + 0.5*std::pow(-mym_BBC3 + ym_BBC3, 2)/std::pow(sigma_ym_BBC3, 2);
            break;
        case 826:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_PMAIP1, 2)) + 0.5*std::pow(-mym_PMAIP1 + ym_PMAIP1, 2)/std::pow(sigma_ym_PMAIP1, 2);
            break;
        case 827:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_BCL2L11, 2)) + 0.5*std::pow(-mym_BCL2L11 + ym_BCL2L11, 2)/std::pow(sigma_ym_BCL2L11, 2);
            break;
        case 828:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_EGF, 2)) + 0.5*std::pow(-mym_EGF + ym_EGF, 2)/std::pow(sigma_ym_EGF, 2);
            break;
        case 829:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_NRG1, 2)) + 0.5*std::pow(-mym_NRG1 + ym_NRG1, 2)/std::pow(sigma_ym_NRG1, 2);
            break;
        case 830:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_EGFR, 2)) + 0.5*std::pow(-mym_EGFR + ym_EGFR, 2)/std::pow(sigma_ym_EGFR, 2);
            break;
        case 831:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_ERBB2, 2)) + 0.5*std::pow(-mym_ERBB2 + ym_ERBB2, 2)/std::pow(sigma_ym_ERBB2, 2);
            break;
        case 832:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_ERBB3, 2)) + 0.5*std::pow(-mym_ERBB3 + ym_ERBB3, 2)/std::pow(sigma_ym_ERBB3, 2);
            break;
        case 833:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_ERBB4, 2)) + 0.5*std::pow(-mym_ERBB4 + ym_ERBB4, 2)/std::pow(sigma_ym_ERBB4, 2);
            break;
        case 834:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_EGFRvIII, 2)) + 0.5*std::pow(-mym_EGFRvIII + ym_EGFRvIII, 2)/std::pow(sigma_ym_EGFRvIII, 2);
            break;
        case 835:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_MET, 2)) + 0.5*std::pow(-mym_MET + ym_MET, 2)/std::pow(sigma_ym_MET, 2);
            break;
        case 836:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_HGF, 2)) + 0.5*std::pow(-mym_HGF + ym_HGF, 2)/std::pow(sigma_ym_HGF, 2);
            break;
        case 837:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_PDGFRA, 2)) + 0.5*std::pow(-mym_PDGFRA + ym_PDGFRA, 2)/std::pow(sigma_ym_PDGFRA, 2);
            break;
        case 838:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_PDGFRB, 2)) + 0.5*std::pow(-mym_PDGFRB + ym_PDGFRB, 2)/std::pow(sigma_ym_PDGFRB, 2);
            break;
        case 839:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_PDGFB, 2)) + 0.5*std::pow(-mym_PDGFB + ym_PDGFB, 2)/std::pow(sigma_ym_PDGFB, 2);
            break;
        case 840:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_SPRY2, 2)) + 0.5*std::pow(-mym_SPRY2 + ym_SPRY2, 2)/std::pow(sigma_ym_SPRY2, 2);
            break;
        case 841:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_CBL, 2)) + 0.5*std::pow(-mym_CBL + ym_CBL, 2)/std::pow(sigma_ym_CBL, 2);
            break;
        case 842:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_GRB2, 2)) + 0.5*std::pow(-mym_GRB2 + ym_GRB2, 2)/std::pow(sigma_ym_GRB2, 2);
            break;
        case 843:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_PLCG1, 2)) + 0.5*std::pow(-mym_PLCG1 + ym_PLCG1, 2)/std::pow(sigma_ym_PLCG1, 2);
            break;
        case 844:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_PLCG2, 2)) + 0.5*std::pow(-mym_PLCG2 + ym_PLCG2, 2)/std::pow(sigma_ym_PLCG2, 2);
            break;
        case 845:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_PIK3CA, 2)) + 0.5*std::pow(-mym_PIK3CA + ym_PIK3CA, 2)/std::pow(sigma_ym_PIK3CA, 2);
            break;
        case 846:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_PIK3CB, 2)) + 0.5*std::pow(-mym_PIK3CB + ym_PIK3CB, 2)/std::pow(sigma_ym_PIK3CB, 2);
            break;
        case 847:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_PIK3CG, 2)) + 0.5*std::pow(-mym_PIK3CG + ym_PIK3CG, 2)/std::pow(sigma_ym_PIK3CG, 2);
            break;
        case 848:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_PIK3CD, 2)) + 0.5*std::pow(-mym_PIK3CD + ym_PIK3CD, 2)/std::pow(sigma_ym_PIK3CD, 2);
            break;
        case 849:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_PIK3R1, 2)) + 0.5*std::pow(-mym_PIK3R1 + ym_PIK3R1, 2)/std::pow(sigma_ym_PIK3R1, 2);
            break;
        case 850:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_PIK3R2, 2)) + 0.5*std::pow(-mym_PIK3R2 + ym_PIK3R2, 2)/std::pow(sigma_ym_PIK3R2, 2);
            break;
        case 851:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_PIK3R3, 2)) + 0.5*std::pow(-mym_PIK3R3 + ym_PIK3R3, 2)/std::pow(sigma_ym_PIK3R3, 2);
            break;
        case 852:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_PIK3R4, 2)) + 0.5*std::pow(-mym_PIK3R4 + ym_PIK3R4, 2)/std::pow(sigma_ym_PIK3R4, 2);
            break;
        case 853:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_PIK3C2A, 2)) + 0.5*std::pow(-mym_PIK3C2A + ym_PIK3C2A, 2)/std::pow(sigma_ym_PIK3C2A, 2);
            break;
        case 854:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_RASGRP1, 2)) + 0.5*std::pow(-mym_RASGRP1 + ym_RASGRP1, 2)/std::pow(sigma_ym_RASGRP1, 2);
            break;
        case 855:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_RASGRP3, 2)) + 0.5*std::pow(-mym_RASGRP3 + ym_RASGRP3, 2)/std::pow(sigma_ym_RASGRP3, 2);
            break;
        case 856:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_NRAS, 2)) + 0.5*std::pow(-mym_NRAS + ym_NRAS, 2)/std::pow(sigma_ym_NRAS, 2);
            break;
        case 857:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_KRAS, 2)) + 0.5*std::pow(-mym_KRAS + ym_KRAS, 2)/std::pow(sigma_ym_KRAS, 2);
            break;
        case 858:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_HRAS, 2)) + 0.5*std::pow(-mym_HRAS + ym_HRAS, 2)/std::pow(sigma_ym_HRAS, 2);
            break;
        case 859:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_NF1, 2)) + 0.5*std::pow(-mym_NF1 + ym_NF1, 2)/std::pow(sigma_ym_NF1, 2);
            break;
        case 860:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_RAF1, 2)) + 0.5*std::pow(-mym_RAF1 + ym_RAF1, 2)/std::pow(sigma_ym_RAF1, 2);
            break;
        case 861:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_BRAF, 2)) + 0.5*std::pow(-mym_BRAF + ym_BRAF, 2)/std::pow(sigma_ym_BRAF, 2);
            break;
        case 862:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_MAP2K1, 2)) + 0.5*std::pow(-mym_MAP2K1 + ym_MAP2K1, 2)/std::pow(sigma_ym_MAP2K1, 2);
            break;
        case 863:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_MAP2K2, 2)) + 0.5*std::pow(-mym_MAP2K2 + ym_MAP2K2, 2)/std::pow(sigma_ym_MAP2K2, 2);
            break;
        case 864:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_DUSP6, 2)) + 0.5*std::pow(-mym_DUSP6 + ym_DUSP6, 2)/std::pow(sigma_ym_DUSP6, 2);
            break;
        case 865:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_RPS6KA1, 2)) + 0.5*std::pow(-mym_RPS6KA1 + ym_RPS6KA1, 2)/std::pow(sigma_ym_RPS6KA1, 2);
            break;
        case 866:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_RPS6KA2, 2)) + 0.5*std::pow(-mym_RPS6KA2 + ym_RPS6KA2, 2)/std::pow(sigma_ym_RPS6KA2, 2);
            break;
        case 867:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_RPS6KA3, 2)) + 0.5*std::pow(-mym_RPS6KA3 + ym_RPS6KA3, 2)/std::pow(sigma_ym_RPS6KA3, 2);
            break;
        case 868:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_RPS6KA4, 2)) + 0.5*std::pow(-mym_RPS6KA4 + ym_RPS6KA4, 2)/std::pow(sigma_ym_RPS6KA4, 2);
            break;
        case 869:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_DUSP1, 2)) + 0.5*std::pow(-mym_DUSP1 + ym_DUSP1, 2)/std::pow(sigma_ym_DUSP1, 2);
            break;
        case 870:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_FOS, 2)) + 0.5*std::pow(-mym_FOS + ym_FOS, 2)/std::pow(sigma_ym_FOS, 2);
            break;
        case 871:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_JUN, 2)) + 0.5*std::pow(-mym_JUN + ym_JUN, 2)/std::pow(sigma_ym_JUN, 2);
            break;
        case 872:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_MYC, 2)) + 0.5*std::pow(-mym_MYC + ym_MYC, 2)/std::pow(sigma_ym_MYC, 2);
            break;
        case 873:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_CTNNB1, 2)) + 0.5*std::pow(-mym_CTNNB1 + ym_CTNNB1, 2)/std::pow(sigma_ym_CTNNB1, 2);
            break;
        case 874:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_PTEN, 2)) + 0.5*std::pow(-mym_PTEN + ym_PTEN, 2)/std::pow(sigma_ym_PTEN, 2);
            break;
        case 875:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_AKT1, 2)) + 0.5*std::pow(-mym_AKT1 + ym_AKT1, 2)/std::pow(sigma_ym_AKT1, 2);
            break;
        case 876:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_AKT2, 2)) + 0.5*std::pow(-mym_AKT2 + ym_AKT2, 2)/std::pow(sigma_ym_AKT2, 2);
            break;
        case 877:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_PDPK1, 2)) + 0.5*std::pow(-mym_PDPK1 + ym_PDPK1, 2)/std::pow(sigma_ym_PDPK1, 2);
            break;
        case 878:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_RICTOR, 2)) + 0.5*std::pow(-mym_RICTOR + ym_RICTOR, 2)/std::pow(sigma_ym_RICTOR, 2);
            break;
        case 879:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_MTOR, 2)) + 0.5*std::pow(-mym_MTOR + ym_MTOR, 2)/std::pow(sigma_ym_MTOR, 2);
            break;
        case 880:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_GSK3B, 2)) + 0.5*std::pow(-mym_GSK3B + ym_GSK3B, 2)/std::pow(sigma_ym_GSK3B, 2);
            break;
        case 881:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_TSC1, 2)) + 0.5*std::pow(-mym_TSC1 + ym_TSC1, 2)/std::pow(sigma_ym_TSC1, 2);
            break;
        case 882:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_TSC2, 2)) + 0.5*std::pow(-mym_TSC2 + ym_TSC2, 2)/std::pow(sigma_ym_TSC2, 2);
            break;
        case 883:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_PRKCA, 2)) + 0.5*std::pow(-mym_PRKCA + ym_PRKCA, 2)/std::pow(sigma_ym_PRKCA, 2);
            break;
        case 884:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_PRKCB, 2)) + 0.5*std::pow(-mym_PRKCB + ym_PRKCB, 2)/std::pow(sigma_ym_PRKCB, 2);
            break;
        case 885:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_PRKCG, 2)) + 0.5*std::pow(-mym_PRKCG + ym_PRKCG, 2)/std::pow(sigma_ym_PRKCG, 2);
            break;
        case 886:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_PRKCD, 2)) + 0.5*std::pow(-mym_PRKCD + ym_PRKCD, 2)/std::pow(sigma_ym_PRKCD, 2);
            break;
        case 887:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_PEBP1, 2)) + 0.5*std::pow(-mym_PEBP1 + ym_PEBP1, 2)/std::pow(sigma_ym_PEBP1, 2);
            break;
        case 888:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_MAPK1, 2)) + 0.5*std::pow(-mym_MAPK1 + ym_MAPK1, 2)/std::pow(sigma_ym_MAPK1, 2);
            break;
        case 889:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_MAPK3, 2)) + 0.5*std::pow(-mym_MAPK3 + ym_MAPK3, 2)/std::pow(sigma_ym_MAPK3, 2);
            break;
        case 890:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_FOXO3, 2)) + 0.5*std::pow(-mym_FOXO3 + ym_FOXO3, 2)/std::pow(sigma_ym_FOXO3, 2);
            break;
        case 891:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_RHEB, 2)) + 0.5*std::pow(-mym_RHEB + ym_RHEB, 2)/std::pow(sigma_ym_RHEB, 2);
            break;
        case 892:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_RPTOR, 2)) + 0.5*std::pow(-mym_RPTOR + ym_RPTOR, 2)/std::pow(sigma_ym_RPTOR, 2);
            break;
        case 893:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_RPS6KB1, 2)) + 0.5*std::pow(-mym_RPS6KB1 + ym_RPS6KB1, 2)/std::pow(sigma_ym_RPS6KB1, 2);
            break;
        case 894:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_RPS6KB2, 2)) + 0.5*std::pow(-mym_RPS6KB2 + ym_RPS6KB2, 2)/std::pow(sigma_ym_RPS6KB2, 2);
            break;
        case 895:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_EIF4EBP1, 2)) + 0.5*std::pow(-mym_EIF4EBP1 + ym_EIF4EBP1, 2)/std::pow(sigma_ym_EIF4EBP1, 2);
            break;
        case 896:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_SOS1, 2)) + 0.5*std::pow(-mym_SOS1 + ym_SOS1, 2)/std::pow(sigma_ym_SOS1, 2);
            break;
        case 897:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_CDKN2A, 2)) + 0.5*std::pow(-mym_CDKN2A + ym_CDKN2A, 2)/std::pow(sigma_ym_CDKN2A, 2);
            break;
        case 898:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_MDM4, 2)) + 0.5*std::pow(-mym_MDM4 + ym_MDM4, 2)/std::pow(sigma_ym_MDM4, 2);
            break;
        case 899:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_FGFR1, 2)) + 0.5*std::pow(-mym_FGFR1 + ym_FGFR1, 2)/std::pow(sigma_ym_FGFR1, 2);
            break;
        case 900:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_FGFR2, 2)) + 0.5*std::pow(-mym_FGFR2 + ym_FGFR2, 2)/std::pow(sigma_ym_FGFR2, 2);
            break;
        case 901:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_FGF1, 2)) + 0.5*std::pow(-mym_FGF1 + ym_FGF1, 2)/std::pow(sigma_ym_FGF1, 2);
            break;
        case 902:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_FGF2, 2)) + 0.5*std::pow(-mym_FGF2 + ym_FGF2, 2)/std::pow(sigma_ym_FGF2, 2);
            break;
        case 903:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_EIF4E, 2)) + 0.5*std::pow(-mym_EIF4E + ym_EIF4E, 2)/std::pow(sigma_ym_EIF4E, 2);
            break;
        case 904:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_IRS1, 2)) + 0.5*std::pow(-mym_IRS1 + ym_IRS1, 2)/std::pow(sigma_ym_IRS1, 2);
            break;
        case 905:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_IRS2, 2)) + 0.5*std::pow(-mym_IRS2 + ym_IRS2, 2)/std::pow(sigma_ym_IRS2, 2);
            break;
        case 906:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_IGF1, 2)) + 0.5*std::pow(-mym_IGF1 + ym_IGF1, 2)/std::pow(sigma_ym_IGF1, 2);
            break;
        case 907:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_IGF2, 2)) + 0.5*std::pow(-mym_IGF2 + ym_IGF2, 2)/std::pow(sigma_ym_IGF2, 2);
            break;
        case 908:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_IGF1R, 2)) + 0.5*std::pow(-mym_IGF1R + ym_IGF1R, 2)/std::pow(sigma_ym_IGF1R, 2);
            break;
        case 909:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_MSH6, 2)) + 0.5*std::pow(-mym_MSH6 + ym_MSH6, 2)/std::pow(sigma_ym_MSH6, 2);
            break;
        case 910:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_BRCA2, 2)) + 0.5*std::pow(-mym_BRCA2 + ym_BRCA2, 2)/std::pow(sigma_ym_BRCA2, 2);
            break;
        case 911:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_MGMT, 2)) + 0.5*std::pow(-mym_MGMT + ym_MGMT, 2)/std::pow(sigma_ym_MGMT, 2);
            break;
        case 912:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_INSR, 2)) + 0.5*std::pow(-mym_INSR + ym_INSR, 2)/std::pow(sigma_ym_INSR, 2);
            break;
        case 913:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ym_INS, 2)) + 0.5*std::pow(-mym_INS + ym_INS, 2)/std::pow(sigma_ym_INS, 2);
            break;
        case 914:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yCytoplasm, 2)) + 0.5*std::pow(-myCytoplasm + yCytoplasm, 2)/std::pow(sigma_yCytoplasm, 2);
            break;
        case 915:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yExtracellular, 2)) + 0.5*std::pow(-myExtracellular + yExtracellular, 2)/std::pow(sigma_yExtracellular, 2);
            break;
        case 916:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yNucleus, 2)) + 0.5*std::pow(-myNucleus + yNucleus, 2)/std::pow(sigma_yNucleus, 2);
            break;
        case 917:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yMitochondrion, 2)) + 0.5*std::pow(-myMitochondrion + yMitochondrion, 2)/std::pow(sigma_yMitochondrion, 2);
            break;
    }
}

} // namespace amici
} // namespace model_SPARCED