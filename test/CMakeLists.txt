find_package(Boost REQUIRED unit_test_framework)

target_sources("testrunner" # Target name
        PRIVATE # Visibility

        # Included sources
        "testExample.cpp"
)

target_link_libraries("testrunner" PRIVATE Boost::headers Boost::unit_test_framework ${DER_LIBARAY_NAME})
