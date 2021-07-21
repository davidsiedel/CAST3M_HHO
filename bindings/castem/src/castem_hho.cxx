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
hho_equ_ss_tri3_get_element_description(
    castem_hho_element_description *const description) {
  try {
      using Cast3MFunc = feta::hybrid::Cas3MFunctions<ElementType::LinearTriangle, feta::hybrid::FieldType::VECTOR_SMALL_STRAIN, 1, 1>;
      Cast3MFunc::getDescription(description);
  } catch (...) {
    return castem_hho_handle_cxx_exception();
  }
  return castem_hho_report_success();
}

MECHHANO_CASTEM_EXPORT castem_hho_status
hho_equ_ss_tri3_get_element_functions(
        castem_hho_element_functions *const elem_funs) {
    try {
        using Cast3MFunc = feta::hybrid::Cas3MFunctions<ElementType::LinearTriangle, feta::hybrid::FieldType::VECTOR_SMALL_STRAIN, 1, 1>;
        Cast3MFunc::getFunctions(elem_funs);
    } catch (...) {
        return castem_hho_handle_cxx_exception();
    }
    return castem_hho_report_success();
}

MECHHANO_CASTEM_EXPORT castem_hho_status
hho_get_generic_functions(
        castem_hho_generic_functions *const gene_funs) {
    try {
        feta::hybrid::getGenericFunctions(gene_funs);
    } catch (...) {
        return castem_hho_handle_cxx_exception();
    }
    return castem_hho_report_success();
}

}  // end of externe "C"
