# C library containing the basis of Cast3M bindings

mechhano_library(${PROJECT_NAME}Cast3MBindings-c
        Cast3MBindings.cxx)
target_include_directories(${PROJECT_NAME}Cast3MBindings-c
   PUBLIC 
   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/bindings/castem/include>
   $<INSTALL_INTERFACE:include>)
target_link_libraries(${PROJECT_NAME}Cast3MBindings-c PUBLIC ${PROJECT_NAME})
if(NOT WIN32)
  target_link_libraries(${PROJECT_NAME}Cast3MBindings-c PRIVATE dl)
endif(NOT WIN32)

# main fortran library

mechhano_library(${PROJECT_NAME}Cast3MBindings SHARED
		castem_hho_fortran_utilities.f95
		castem_hho.f95)
set_target_properties(${PROJECT_NAME}Cast3MBindings
  PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/bindings/castem/modules")
target_include_directories(${PROJECT_NAME}Cast3MBindings
	PRIVATE "${CMAKE_SOURCE_DIR}/include"
	PRIVATE "${CMAKE_SOURCE_DIR}/bindings/castem/include")
target_link_libraries(${PROJECT_NAME}Cast3MBindings
   PUBLIC ${PROJECT_NAME}Cast3MBindings-c)

# library containing a set of predefined elements
 
mechhano_library(${PROJECT_NAME}Cast3MElements
    castem_hho.cxx)
target_link_libraries(${PROJECT_NAME}Cast3MElements
	PUBLIC ${PROJECT_NAME}Cast3MBindings-c)
