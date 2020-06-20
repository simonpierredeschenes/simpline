cmake_minimum_required(VERSION 2.8.3)
project(simpline)

set(CMAKE_CXX_STANDARD 11)

set(EXTERNAL_INCLUDE_DIRS "")
set(EXTERNAL_LIBS "")

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
set(EXTERNAL_INCLUDE_DIRS ${EXTERNAL_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR})

add_library(simpline SHARED simpline/ParametrizedSpline.cpp simpline/ConstantSpeedSpline.cpp)

target_link_libraries(simpline Eigen3::Eigen)

# install target
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")
set(DEF_INSTALL_CMAKE_DIR lib/cmake/simpline CACHE PATH "Installation directory for CMake files")

install(TARGETS simpline DESTINATION ${INSTALL_LIB_DIR})

install(FILES simpline/Simpline.h
  DESTINATION ${INSTALL_INCLUDE_DIR}/simpline
  )

configure_file(simplineConfig.cmake.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/simplineConfig.cmake" @ONLY
  )
install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/simplineConfig.cmake"
  DESTINATION ${DEF_INSTALL_CMAKE_DIR}
  COMPONENT dev
  )

# uninstall target
if (NOT TARGET uninstall)
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)
  add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif ()
