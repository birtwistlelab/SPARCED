#include "amici/model.h"
#include "wrapfunctions.h"
#include "SPARCED_IFNg.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_SPARCED_IFNg::Model_SPARCED_IFNg());
}


} // namespace generic_model

} // namespace amici
