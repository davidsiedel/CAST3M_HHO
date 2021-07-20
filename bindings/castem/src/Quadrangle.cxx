/*!
 * \file   Triangle.cxx
 * \brief
 * \author Thomas Helfer
 * \date   10/03/2021
 */

#include "Mechhano/Cast3M/Config.hxx"
#include "Mechhano/Cast3M/Cast3MBindings.hxx"

//#include "Mechhano/func.hxx"
#include "feta/functions.hxx"

extern "C" {

// public interface for the element
MECHHANO_CASTEM_EXPORT castem_hho_status
qua4_1_1_hho_ss(
        castem_hho_element_description *const description,
castem_hho_element_geometry *const gdescription) {
try {
const int a = 1;
using Cast3MFunc = feta::hybrid::Cas3MFunctions<ElementType::LinearQuadrangle, feta::hybrid::FieldType::VECTOR_SMALL_STRAIN, 1, 1>;
Cast3MFunc::getDescription(description, gdescription);
} catch (...) {
return castem_hho_handle_cxx_exception();
}
return castem_hho_report_success();
}

}  // end of externe "C"
