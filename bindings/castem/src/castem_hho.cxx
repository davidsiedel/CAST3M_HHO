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
hho_1_1_sp_tri3_get_element_description(
    castem_hho_element_description *const description) {
  try {
      using Cast3MFunc = feta::hybrid::Cas3MFunctions<ElementType::LinearTriangle, feta::hybrid::FieldType::VECTOR_SMALL_STRAIN_CAST3M, 1, 1>;
      Cast3MFunc::getDescription(description);
  } catch (...) {
    return castem_hho_handle_cxx_exception();
  }
  return castem_hho_report_success();
}

MECHHANO_CASTEM_EXPORT castem_hho_status
hho_1_1_sp_tri3_get_element_functions(
        castem_hho_element_functions *const elem_funs) {
    try {
        using Cast3MFunc = feta::hybrid::Cas3MFunctions<ElementType::LinearTriangle, feta::hybrid::FieldType::VECTOR_SMALL_STRAIN_CAST3M, 1, 1>;
        Cast3MFunc::getFunctions(elem_funs);
    } catch (...) {
        return castem_hho_handle_cxx_exception();
    }
    return castem_hho_report_success();
}
MECHHANO_CASTEM_EXPORT castem_hho_status
hho_1_1_sp_qua4_get_element_description(
        castem_hho_element_description *const description) {
    try {
        using Cast3MFunc = feta::hybrid::Cas3MFunctions<ElementType::LinearQuadrangle, feta::hybrid::FieldType::VECTOR_SMALL_STRAIN_CAST3M, 1, 1>;
        Cast3MFunc::getDescription(description);
    } catch (...) {
        return castem_hho_handle_cxx_exception();
    }
    return castem_hho_report_success();
}

MECHHANO_CASTEM_EXPORT castem_hho_status
hho_1_1_sp_qua4_get_element_functions(
        castem_hho_element_functions *const elem_funs) {
    try {
        using Cast3MFunc = feta::hybrid::Cas3MFunctions<ElementType::LinearQuadrangle, feta::hybrid::FieldType::VECTOR_SMALL_STRAIN_CAST3M, 1, 1>;
        Cast3MFunc::getFunctions(elem_funs);
    } catch (...) {
        return castem_hho_handle_cxx_exception();
    }
    return castem_hho_report_success();
}

MECHHANO_CASTEM_EXPORT castem_hho_status
hho_1_1_sp_poly05_get_element_description(
        castem_hho_element_description *const description) {
    try {
        using Cast3MFunc = feta::hybrid::Cas3MFunctions<ElementType::LinearPentagon, feta::hybrid::FieldType::VECTOR_SMALL_STRAIN_CAST3M, 1, 1>;
        Cast3MFunc::getDescription(description);
    } catch (...) {
        return castem_hho_handle_cxx_exception();
    }
    return castem_hho_report_success();
}

MECHHANO_CASTEM_EXPORT castem_hho_status
hho_1_1_sp_poly05_get_element_functions(
        castem_hho_element_functions *const elem_funs) {
    try {
        using Cast3MFunc = feta::hybrid::Cas3MFunctions<ElementType::LinearPentagon, feta::hybrid::FieldType::VECTOR_SMALL_STRAIN_CAST3M, 1, 1>;
        Cast3MFunc::getFunctions(elem_funs);
    } catch (...) {
        return castem_hho_handle_cxx_exception();
    }
    return castem_hho_report_success();
}

MECHHANO_CASTEM_EXPORT castem_hho_status
hho_1_1_sp_poly06_get_element_description(
    castem_hho_element_description *const description) {
  try {
    using Cast3MFunc = feta::hybrid::Cas3MFunctions<ElementType::LinearHexagon, feta::hybrid::FieldType::VECTOR_SMALL_STRAIN_CAST3M, 1, 1>;
    Cast3MFunc::getDescription(description);
  } catch (...) {
    return castem_hho_handle_cxx_exception();
  }
  return castem_hho_report_success();
}
MECHHANO_CASTEM_EXPORT castem_hho_status
hho_1_1_sp_poly06_get_element_functions(
    castem_hho_element_functions *const elem_funs) {
  try {
    using Cast3MFunc = feta::hybrid::Cas3MFunctions<ElementType::LinearHexagon, feta::hybrid::FieldType::VECTOR_SMALL_STRAIN_CAST3M, 1, 1>;
    Cast3MFunc::getFunctions(elem_funs);
  } catch (...) {
    return castem_hho_handle_cxx_exception();
  }
  return castem_hho_report_success();
}

MECHHANO_CASTEM_EXPORT castem_hho_status
hho_1_1_sp_poly07_get_element_description(
    castem_hho_element_description *const description) {
  try {
    using Cast3MFunc = feta::hybrid::Cas3MFunctions<ElementType::LinearHeptagon, feta::hybrid::FieldType::VECTOR_SMALL_STRAIN_CAST3M, 1, 1>;
    Cast3MFunc::getDescription(description);
  } catch (...) {
    return castem_hho_handle_cxx_exception();
  }
  return castem_hho_report_success();
}
MECHHANO_CASTEM_EXPORT castem_hho_status
hho_1_1_sp_poly07_get_element_functions(
    castem_hho_element_functions *const elem_funs) {
  try {
    using Cast3MFunc = feta::hybrid::Cas3MFunctions<ElementType::LinearHeptagon, feta::hybrid::FieldType::VECTOR_SMALL_STRAIN_CAST3M, 1, 1>;
    Cast3MFunc::getFunctions(elem_funs);
  } catch (...) {
    return castem_hho_handle_cxx_exception();
  }
  return castem_hho_report_success();
}

MECHHANO_CASTEM_EXPORT castem_hho_status
hho_1_1_sp_poly08_get_element_description(
    castem_hho_element_description *const description) {
  try {
    using Cast3MFunc = feta::hybrid::Cas3MFunctions<ElementType::LinearOctagon, feta::hybrid::FieldType::VECTOR_SMALL_STRAIN_CAST3M, 1, 1>;
    Cast3MFunc::getDescription(description);
  } catch (...) {
    return castem_hho_handle_cxx_exception();
  }
  return castem_hho_report_success();
}
MECHHANO_CASTEM_EXPORT castem_hho_status
hho_1_1_sp_poly08_get_element_functions(
    castem_hho_element_functions *const elem_funs) {
  try {
    using Cast3MFunc = feta::hybrid::Cas3MFunctions<ElementType::LinearOctagon, feta::hybrid::FieldType::VECTOR_SMALL_STRAIN_CAST3M, 1, 1>;
    Cast3MFunc::getFunctions(elem_funs);
  } catch (...) {
    return castem_hho_handle_cxx_exception();
  }
  return castem_hho_report_success();
}

MECHHANO_CASTEM_EXPORT castem_hho_status
hho_1_1_sp_poly09_get_element_description(
    castem_hho_element_description *const description) {
  try {
    using Cast3MFunc = feta::hybrid::Cas3MFunctions<ElementType::LinearEnneagon, feta::hybrid::FieldType::VECTOR_SMALL_STRAIN_CAST3M, 1, 1>;
    Cast3MFunc::getDescription(description);
  } catch (...) {
    return castem_hho_handle_cxx_exception();
  }
  return castem_hho_report_success();
}
MECHHANO_CASTEM_EXPORT castem_hho_status
hho_1_1_sp_poly09_get_element_functions(
    castem_hho_element_functions *const elem_funs) {
  try {
    using Cast3MFunc = feta::hybrid::Cas3MFunctions<ElementType::LinearEnneagon, feta::hybrid::FieldType::VECTOR_SMALL_STRAIN_CAST3M, 1, 1>;
    Cast3MFunc::getFunctions(elem_funs);
  } catch (...) {
    return castem_hho_handle_cxx_exception();
  }
  return castem_hho_report_success();
}

MECHHANO_CASTEM_EXPORT castem_hho_status
hho_1_1_sp_poly10_get_element_description(
    castem_hho_element_description *const description) {
  try {
    using Cast3MFunc = feta::hybrid::Cas3MFunctions<ElementType::LinearDecagon, feta::hybrid::FieldType::VECTOR_SMALL_STRAIN_CAST3M, 1, 1>;
    Cast3MFunc::getDescription(description);
  } catch (...) {
    return castem_hho_handle_cxx_exception();
  }
  return castem_hho_report_success();
}
MECHHANO_CASTEM_EXPORT castem_hho_status
hho_1_1_sp_poly10_get_element_functions(
    castem_hho_element_functions *const elem_funs) {
  try {
    using Cast3MFunc = feta::hybrid::Cas3MFunctions<ElementType::LinearDecagon, feta::hybrid::FieldType::VECTOR_SMALL_STRAIN_CAST3M, 1, 1>;
    Cast3MFunc::getFunctions(elem_funs);
  } catch (...) {
    return castem_hho_handle_cxx_exception();
  }
  return castem_hho_report_success();
}

MECHHANO_CASTEM_EXPORT castem_hho_status
hho_1_1_sp_poly11_get_element_description(
    castem_hho_element_description *const description) {
  try {
    using Cast3MFunc = feta::hybrid::Cas3MFunctions<ElementType::LinearHendecagon, feta::hybrid::FieldType::VECTOR_SMALL_STRAIN_CAST3M, 1, 1>;
    Cast3MFunc::getDescription(description);
  } catch (...) {
    return castem_hho_handle_cxx_exception();
  }
  return castem_hho_report_success();
}
MECHHANO_CASTEM_EXPORT castem_hho_status
hho_1_1_sp_poly11_get_element_functions(
    castem_hho_element_functions *const elem_funs) {
  try {
    using Cast3MFunc = feta::hybrid::Cas3MFunctions<ElementType::LinearHendecagon, feta::hybrid::FieldType::VECTOR_SMALL_STRAIN_CAST3M, 1, 1>;
    Cast3MFunc::getFunctions(elem_funs);
  } catch (...) {
    return castem_hho_handle_cxx_exception();
  }
  return castem_hho_report_success();
}

MECHHANO_CASTEM_EXPORT castem_hho_status
hho_1_1_sp_poly12_get_element_description(
    castem_hho_element_description *const description) {
  try {
    using Cast3MFunc = feta::hybrid::Cas3MFunctions<ElementType::LinearDodecagon, feta::hybrid::FieldType::VECTOR_SMALL_STRAIN_CAST3M, 1, 1>;
    Cast3MFunc::getDescription(description);
  } catch (...) {
    return castem_hho_handle_cxx_exception();
  }
  return castem_hho_report_success();
}
MECHHANO_CASTEM_EXPORT castem_hho_status
hho_1_1_sp_poly12_get_element_functions(
    castem_hho_element_functions *const elem_funs) {
  try {
    using Cast3MFunc = feta::hybrid::Cas3MFunctions<ElementType::LinearDodecagon, feta::hybrid::FieldType::VECTOR_SMALL_STRAIN_CAST3M, 1, 1>;
    Cast3MFunc::getFunctions(elem_funs);
  } catch (...) {
    return castem_hho_handle_cxx_exception();
  }
  return castem_hho_report_success();
}

MECHHANO_CASTEM_EXPORT castem_hho_status
hho_1_1_sp_poly13_get_element_description(
    castem_hho_element_description *const description) {
  try {
    using Cast3MFunc = feta::hybrid::Cas3MFunctions<ElementType::LinearTridecagon, feta::hybrid::FieldType::VECTOR_SMALL_STRAIN_CAST3M, 1, 1>;
    Cast3MFunc::getDescription(description);
  } catch (...) {
    return castem_hho_handle_cxx_exception();
  }
  return castem_hho_report_success();
}
MECHHANO_CASTEM_EXPORT castem_hho_status
hho_1_1_sp_poly13_get_element_functions(
    castem_hho_element_functions *const elem_funs) {
  try {
    using Cast3MFunc = feta::hybrid::Cas3MFunctions<ElementType::LinearTridecagon, feta::hybrid::FieldType::VECTOR_SMALL_STRAIN_CAST3M, 1, 1>;
    Cast3MFunc::getFunctions(elem_funs);
  } catch (...) {
    return castem_hho_handle_cxx_exception();
  }
  return castem_hho_report_success();
}

MECHHANO_CASTEM_EXPORT castem_hho_status
hho_1_1_sp_poly14_get_element_description(
    castem_hho_element_description *const description) {
  try {
    using Cast3MFunc = feta::hybrid::Cas3MFunctions<ElementType::LinearTetradecagon, feta::hybrid::FieldType::VECTOR_SMALL_STRAIN_CAST3M, 1, 1>;
    Cast3MFunc::getDescription(description);
  } catch (...) {
    return castem_hho_handle_cxx_exception();
  }
  return castem_hho_report_success();
}
MECHHANO_CASTEM_EXPORT castem_hho_status
hho_1_1_sp_poly14_get_element_functions(
    castem_hho_element_functions *const elem_funs) {
  try {
    using Cast3MFunc = feta::hybrid::Cas3MFunctions<ElementType::LinearTetradecagon, feta::hybrid::FieldType::VECTOR_SMALL_STRAIN_CAST3M, 1, 1>;
    Cast3MFunc::getFunctions(elem_funs);
  } catch (...) {
    return castem_hho_handle_cxx_exception();
  }
  return castem_hho_report_success();
}

MECHHANO_CASTEM_EXPORT castem_hho_status
hho_1_1_sp_poly15_get_element_description(
    castem_hho_element_description *const description) {
  try {
    using Cast3MFunc = feta::hybrid::Cas3MFunctions<ElementType::LinearPentadecagon, feta::hybrid::FieldType::VECTOR_SMALL_STRAIN_CAST3M, 1, 1>;
    Cast3MFunc::getDescription(description);
  } catch (...) {
    return castem_hho_handle_cxx_exception();
  }
  return castem_hho_report_success();
}
MECHHANO_CASTEM_EXPORT castem_hho_status
hho_1_1_sp_poly15_get_element_functions(
    castem_hho_element_functions *const elem_funs) {
  try {
    using Cast3MFunc = feta::hybrid::Cas3MFunctions<ElementType::LinearPentadecagon, feta::hybrid::FieldType::VECTOR_SMALL_STRAIN_CAST3M, 1, 1>;
    Cast3MFunc::getFunctions(elem_funs);
  } catch (...) {
    return castem_hho_handle_cxx_exception();
  }
  return castem_hho_report_success();
}

MECHHANO_CASTEM_EXPORT castem_hho_status
hho_1_1_sp_poly16_get_element_description(
    castem_hho_element_description *const description) {
  try {
    using Cast3MFunc = feta::hybrid::Cas3MFunctions<ElementType::LinearHexadecagon, feta::hybrid::FieldType::VECTOR_SMALL_STRAIN_CAST3M, 1, 1>;
    Cast3MFunc::getDescription(description);
  } catch (...) {
    return castem_hho_handle_cxx_exception();
  }
  return castem_hho_report_success();
}
MECHHANO_CASTEM_EXPORT castem_hho_status
hho_1_1_sp_poly16_get_element_functions(
    castem_hho_element_functions *const elem_funs) {
  try {
    using Cast3MFunc = feta::hybrid::Cas3MFunctions<ElementType::LinearHexadecagon, feta::hybrid::FieldType::VECTOR_SMALL_STRAIN_CAST3M, 1, 1>;
    Cast3MFunc::getFunctions(elem_funs);
  } catch (...) {
    return castem_hho_handle_cxx_exception();
  }
  return castem_hho_report_success();
}

MECHHANO_CASTEM_EXPORT castem_hho_status
hho_1_1_sp_poly17_get_element_description(
    castem_hho_element_description *const description) {
  try {
    using Cast3MFunc = feta::hybrid::Cas3MFunctions<ElementType::LinearHeptadecagon, feta::hybrid::FieldType::VECTOR_SMALL_STRAIN_CAST3M, 1, 1>;
    Cast3MFunc::getDescription(description);
  } catch (...) {
    return castem_hho_handle_cxx_exception();
  }
  return castem_hho_report_success();
}
MECHHANO_CASTEM_EXPORT castem_hho_status
hho_1_1_sp_poly17_get_element_functions(
    castem_hho_element_functions *const elem_funs) {
  try {
    using Cast3MFunc = feta::hybrid::Cas3MFunctions<ElementType::LinearHeptadecagon, feta::hybrid::FieldType::VECTOR_SMALL_STRAIN_CAST3M, 1, 1>;
    Cast3MFunc::getFunctions(elem_funs);
  } catch (...) {
    return castem_hho_handle_cxx_exception();
  }
  return castem_hho_report_success();
}

MECHHANO_CASTEM_EXPORT castem_hho_status
hho_1_1_sp_poly18_get_element_description(
    castem_hho_element_description *const description) {
  try {
    using Cast3MFunc = feta::hybrid::Cas3MFunctions<ElementType::LinearOctadecagon, feta::hybrid::FieldType::VECTOR_SMALL_STRAIN_CAST3M, 1, 1>;
    Cast3MFunc::getDescription(description);
  } catch (...) {
    return castem_hho_handle_cxx_exception();
  }
  return castem_hho_report_success();
}
MECHHANO_CASTEM_EXPORT castem_hho_status
hho_1_1_sp_poly18_get_element_functions(
    castem_hho_element_functions *const elem_funs) {
  try {
    using Cast3MFunc = feta::hybrid::Cas3MFunctions<ElementType::LinearOctadecagon, feta::hybrid::FieldType::VECTOR_SMALL_STRAIN_CAST3M, 1, 1>;
    Cast3MFunc::getFunctions(elem_funs);
  } catch (...) {
    return castem_hho_handle_cxx_exception();
  }
  return castem_hho_report_success();
}

MECHHANO_CASTEM_EXPORT castem_hho_status
hho_1_1_sp_poly19_get_element_description(
    castem_hho_element_description *const description) {
  try {
    using Cast3MFunc = feta::hybrid::Cas3MFunctions<ElementType::LinearEnneadecagon, feta::hybrid::FieldType::VECTOR_SMALL_STRAIN_CAST3M, 1, 1>;
    Cast3MFunc::getDescription(description);
  } catch (...) {
    return castem_hho_handle_cxx_exception();
  }
  return castem_hho_report_success();
}
MECHHANO_CASTEM_EXPORT castem_hho_status
hho_1_1_sp_poly19_get_element_functions(
    castem_hho_element_functions *const elem_funs) {
  try {
    using Cast3MFunc = feta::hybrid::Cas3MFunctions<ElementType::LinearEnneadecagon, feta::hybrid::FieldType::VECTOR_SMALL_STRAIN_CAST3M, 1, 1>;
    Cast3MFunc::getFunctions(elem_funs);
  } catch (...) {
    return castem_hho_handle_cxx_exception();
  }
  return castem_hho_report_success();
}

MECHHANO_CASTEM_EXPORT castem_hho_status
hho_1_1_sp_poly20_get_element_description(
    castem_hho_element_description *const description) {
  try {
    using Cast3MFunc = feta::hybrid::Cas3MFunctions<ElementType::LinearIcosagon, feta::hybrid::FieldType::VECTOR_SMALL_STRAIN_CAST3M, 1, 1>;
    Cast3MFunc::getDescription(description);
  } catch (...) {
    return castem_hho_handle_cxx_exception();
  }
  return castem_hho_report_success();
}
MECHHANO_CASTEM_EXPORT castem_hho_status
hho_1_1_sp_poly20_get_element_functions(
    castem_hho_element_functions *const elem_funs) {
  try {
    using Cast3MFunc = feta::hybrid::Cas3MFunctions<ElementType::LinearIcosagon, feta::hybrid::FieldType::VECTOR_SMALL_STRAIN_CAST3M, 1, 1>;
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
