macro(castem_hho_fortran_module file)
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${file}"
    DESTINATION "include")
endmacro(castem_hho_fortran_module)

castem_hho_fortran_module(castem_hho_fortran_utilities.mod)
castem_hho_fortran_module(castem_hho.mod)
