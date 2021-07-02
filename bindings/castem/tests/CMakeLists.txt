if((CMAKE_HOST_WIN32) AND (NOT MSYS))
	#  set_property(TEST Cast3MTriangleTests
	#    PROPERTY DEPENDS BehaviourTest
	#    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;${MGIS_PATH_STRING}")
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_executable(Cast3MTriangleTests-f Cast3MTriangleTests-f.f95)
target_include_directories(Cast3MTriangleTests-f
  PRIVATE "${CMAKE_BINARY_DIR}/bindings/castem/modules")
target_link_libraries(Cast3MTriangleTests-f ${PROJECT_NAME}Cast3MBindings)
add_test(
  NAME Cast3MTriangleTests_f
  COMMAND Cast3MTriangleTests-f $<TARGET_FILE:${PROJECT_NAME}Cast3MElements>)
set_property(TEST Cast3MTriangleTests_f
    PROPERTY DEPENDS ${PROJECT_NAME}Cast3MElements>)
