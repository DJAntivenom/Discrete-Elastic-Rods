#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-function"
#include <polyscope/polyscope.h>
#pragma GCC diagnostic pop

#include <PolyscopeCallback.hpp>

int main()
{
    overrideCoutBuffer();

    polyscope::options::autocenterStructures = false;
    polyscope::options::autoscaleStructures = false;
    polyscope::options::automaticallyComputeSceneExtents = false;
    polyscope::view::farClipRatio = 50;
    polyscope::state::lengthScale = 1.;
    polyscope::view::setUpDir(polyscope::view::UpDir::YUp);
    polyscope::view::windowWidth = 1920;
    polyscope::view::windowHeight = 1080;
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::TileReflection;
    polyscope::options::groundPlaneHeightFactor = 0.;
    polyscope::options::shadowDarkness = 0.4;
    polyscope::options::buildGui = false;
    polyscope::options::ssaaFactor = 2;
    polyscope::init();

    polyscope::state::userCallback = polyscopeCallback;

    polyscope::show();

    return 0;
}
