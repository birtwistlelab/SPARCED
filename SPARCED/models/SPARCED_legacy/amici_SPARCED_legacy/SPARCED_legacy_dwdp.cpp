#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <array>

#include "tcl.h"
#include "p.h"
#include "k.h"
#include "w.h"
#include "h.h"
#include "x.h"
#include "dtcldp.h"
#include "dwdp.h"

namespace amici {
namespace model_SPARCED_legacy {

void dwdp_SPARCED_legacy(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp){
}

} // namespace amici
} // namespace model_SPARCED_legacy