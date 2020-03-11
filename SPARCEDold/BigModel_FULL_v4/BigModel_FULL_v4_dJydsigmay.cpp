#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "y.h"
#include "my.h"
#include "p.h"
#include "k.h"
#include "sigmay.h"

void dJydsigmay_BigModel_FULL_v4(realtype *dJydsigmay, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigmay[0] = 1.0/sigmaactp53 - 1.0*pow(actp53 - mactp53, 2)/pow(sigmaactp53, 3);
            break;
        case 1:
            dJydsigmay[1] = 1.0/sigmaphERK - 1.0*pow(mphERK - phERK, 2)/pow(sigmaphERK, 3);
            break;
        case 2:
            dJydsigmay[2] = 1.0/sigmaegfLR - 1.0*pow(egfLR - megfLR, 2)/pow(sigmaegfLR, 3);
            break;
    }
}