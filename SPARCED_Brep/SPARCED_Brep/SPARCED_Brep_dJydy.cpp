#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "y.h"
#include "my.h"
#include "p.h"
#include "k.h"
#include "sigmay.h"

void dJydy_SPARCED_Brep(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = 0.5*(2*actp53 - 2*mactp53)/pow(sigmaactp53, 2);
            break;
        case 1:
            dJydy[0] = 0.5*(-2*mphERK + 2*phERK)/pow(sigmaphERK, 2);
            break;
        case 2:
            dJydy[0] = 0.5*(2*egfLR - 2*megfLR)/pow(sigmaegfLR, 2);
            break;
    }
}