#include "amici/model.h"
#include "wrapfunctions.h"
#include "SPARCED_au565.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_SPARCED_au565::Model_SPARCED_au565());
}


} // namespace generic_model

} // namespace amici
