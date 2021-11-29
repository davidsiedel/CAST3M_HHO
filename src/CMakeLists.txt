include(${CMAKE_SOURCE_DIR}/cmake/modules/mechhano.cmake)

mechhano_library(${PROJECT_NAME}
  mechhano.cxx)
target_include_directories(${PROJECT_NAME}
   PUBLIC 
   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
   $<INSTALL_INTERFACE:include>)
target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen)

add_executable(${PROJECT_NAME}_exe
        mechhano.cxx
        )
target_link_libraries(${PROJECT_NAME}_exe PUBLIC ${PROJECT_NAME})
