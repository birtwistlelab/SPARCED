#include "amici/model.h"
#include "wrapfunctions.h"
#include "SPARCED_legacy.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_SPARCED_legacy::Model_SPARCED_legacy());
}


} // namespace generic_model

} // namespace amici
