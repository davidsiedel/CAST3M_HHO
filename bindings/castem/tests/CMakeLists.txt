if((CMAKE_HOST_WIN32) AND (NOT MSYS))
	#  set_property(TEST Cast3MTriangleTests
	#    PROPERTY DEPENDS BehaviourTest
	#    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;${MGIS_PATH_STRING}")
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_executable(test_tri3-f test_triangle.f95)
target_include_directories(test_tri3-f
		PRIVATE "${CMAKE_BINARY_DIR}/bindings/castem/modules")
target_link_libraries(test_tri3-f ${PROJECT_NAME}Cast3MBindings)
add_test(
		NAME test_tri3_f
		COMMAND test_tri3-f $<TARGET_FILE:${PROJECT_NAME}Cast3MElements>)
set_property(TEST test_tri3_f
		PROPERTY DEPENDS ${PROJECT_NAME}Cast3MElements>)

add_executable(test_qua4-f test_quadrangle.f95)
target_include_directories(test_qua4-f
		PRIVATE "${CMAKE_BINARY_DIR}/bindings/castem/modules")
target_link_libraries(test_qua4-f ${PROJECT_NAME}Cast3MBindings)
add_test(
		NAME test_qua4_f
		COMMAND test_qua4-f $<TARGET_FILE:${PROJECT_NAME}Cast3MElements>)
set_property(TEST test_qua4_f
		PROPERTY DEPENDS ${PROJECT_NAME}Cast3MElements>)

add_executable(test_poly05-f test_pentagon.f95)
target_include_directories(test_poly05-f
		PRIVATE "${CMAKE_BINARY_DIR}/bindings/castem/modules")
target_link_libraries(test_poly05-f ${PROJECT_NAME}Cast3MBindings)
add_test(
		NAME test_poly05_f
		COMMAND test_poly05-f $<TARGET_FILE:${PROJECT_NAME}Cast3MElements>)
set_property(TEST test_poly05_f
		PROPERTY DEPENDS ${PROJECT_NAME}Cast3MElements>)

add_executable(test_poly20-f test_icosagon.f95)
target_include_directories(test_poly20-f
		PRIVATE "${CMAKE_BINARY_DIR}/bindings/castem/modules")
target_link_libraries(test_poly20-f ${PROJECT_NAME}Cast3MBindings)
add_test(
		NAME test_poly20_f
		COMMAND test_poly20-f $<TARGET_FILE:${PROJECT_NAME}Cast3MElements>)
set_property(TEST test_poly20_f
		PROPERTY DEPENDS ${PROJECT_NAME}Cast3MElements>)
