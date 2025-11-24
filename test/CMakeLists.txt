find_package(Boost QUIET COMPONENTS unit_test_framework)

if (NOT Boost_FOUND)
  message(STATUS "Tests will not be built because boost was not found")
  return()
elseif(NOT TARGET Boost::unit_test_framework)
  message(STATUS "Tests will not be built because boost is installed without the Boost::unit_test_framework target")
  return()
endif()

# Define test executable
add_executable("testrunner")

# make sure built binary is placed in root directory of build-tree
set_target_properties("testrunner" PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")

target_sources("testrunner" # Target name
        PRIVATE # Visibility

        # Included sources
        "testExample.cpp"
)

target_link_libraries("testrunner" PRIVATE Boost::headers Boost::unit_test_framework ${DER_LIBRARY_NAME})
