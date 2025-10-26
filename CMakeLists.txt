# project definition
cmake_minimum_required(VERSION 3.31)
project("Discrete Elastic Rods" VERSION 0.1.0 LANGUAGES C CXX)

# include helpers
include("cmake/AddCompileOption.cmake")

# add static library. This is the core of the project
set(DER_TARGET_NAME "discrete_elastic_rods")
add_executable(${DER_TARGET_NAME})

# set compile features
target_compile_features(${DER_TARGET_NAME} PUBLIC c_std_23)
target_compile_features(${DER_TARGET_NAME} PUBLIC cxx_std_17)
set_target_properties(${DER_TARGET_NAME} PROPERTIES C_EXTENSIONS OFF CXX_EXTENSIONS OFF)
enable_supported_compiler_flag(${DER_TARGET_NAME} FLAGS "-Wall" "-Wextra")

# setup project
add_subdirectory("src")

# Force colored output for the ninja generator
if(CMAKE_GENERATOR STREQUAL "Ninja")
  enable_supported_compiler_flag(${DER_TARGET_NAME}
                                 PUBLIC
                                 FLAGS "-fcolor-diagnostics" "-fdiagnostics-color=always"
  )
  enable_supported_compiler_flag(${DER_TARGET_NAME}
                                 PUBLIC CHECK_C_FLAGS
                                 FLAGS "-fcolor-diagnostics" "-fdiagnostics-color=always"
  )
endif()
