#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "y.h"
#include "my.h"
#include "p.h"
#include "k.h"
#include "sigmay.h"

void dJydy_BigModel_FULL_v4(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = 1.0*(actp53 - mactp53)/pow(sigmaactp53, 2);
            break;
        case 1:
            dJydy[0] = 1.0*(-mphERK + phERK)/pow(sigmaphERK, 2);
            break;
        case 2:
            dJydy[0] = 1.0*(egfLR - megfLR)/pow(sigmaegfLR, 2);
            break;
    }
}