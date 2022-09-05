#include "amici/model.h"
#include "wrapfunctions.h"
#include "SPARCEDo4a_v1.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_SPARCEDo4a_v1::Model_SPARCEDo4a_v1());
}


} // namespace generic_model

} // namespace amici
