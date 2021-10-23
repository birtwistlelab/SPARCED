#include "amici/model.h"
#include "wrapfunctions.h"
#include "SPARCED_u87i.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_SPARCED_u87i::Model_SPARCED_u87i());
}


} // namespace generic_model

} // namespace amici
