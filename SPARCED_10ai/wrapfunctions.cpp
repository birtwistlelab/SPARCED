#include "amici/model.h"
#include "wrapfunctions.h"
#include "SPARCED_10ai.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_SPARCED_10ai::Model_SPARCED_10ai());
}


} // namespace generic_model

} // namespace amici
