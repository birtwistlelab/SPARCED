#include "amici/model.h"
#include "wrapfunctions.h"
#include "SPARCED.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_SPARCED::Model_SPARCED());
}


} // namespace generic_model

} // namespace amici
