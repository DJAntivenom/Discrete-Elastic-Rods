# find required libraries
include(FetchContent)

FetchContent_Declare(
  polyscope
  GIT_REPOSITORY https://github.com/nmwsharp/polyscope
  GIT_TAG        b1f5ae82488b82d91a61d7d8c81882989216056d # sha of v2.5.0, latest stable release as of 26.10.25
)
FetchContent_Declare(
  eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen
  GIT_TAG        549bf8c75b6aae071cde2f28aa48f16ee3ae60b0 # sha of v5.0.0, latest stable release as of 26.10.25
)

FetchContent_MakeAvailable(polyscope eigen)

# setup executable components
## add include directories
target_include_directories(${DER_TARGET_NAME}
  PUBLIC
  "${PROJECT_SOURCE_DIR}/include"
)

## add dependencies
target_link_libraries(${DER_TARGET_NAME}
  PRIVATE
  polyscope
  Eigen3::Eigen
)

## add sources
target_sources(${DER_TARGET_NAME}
  PRIVATE
  "main.cpp"

  "DiscreteElasticRod.cpp"
  "PolyscopeCallback.cpp"
)
