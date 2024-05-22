#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <array>

#include "p.h"
#include "k.h"
#include "w.h"
#include "h.h"
#include "x.h"
#include "dtcldp.h"

namespace amici {
namespace model_SPARCED_tutorial {

void dydp_SPARCED_tutorial(realtype *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *dtcldp){
}

} // namespace amici
} // namespace model_SPARCED_tutorial