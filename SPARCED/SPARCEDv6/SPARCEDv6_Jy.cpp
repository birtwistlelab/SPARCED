#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "y.h"
#include "my.h"
#include "p.h"
#include "k.h"
#include "sigmay.h"

void Jy_SPARCEDv6(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmaactp53, 2)) + 0.5*pow(actp53 - mactp53, 2)/pow(sigmaactp53, 2);
            break;
        case 1:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmaphERK, 2)) + 0.5*pow(-mphERK + phERK, 2)/pow(sigmaphERK, 2);
            break;
        case 2:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmaegfLR, 2)) + 0.5*pow(egfLR - megfLR, 2)/pow(sigmaegfLR, 2);
            break;
    }
}