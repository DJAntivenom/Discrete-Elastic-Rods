/* Get rid of annoying unused ... warnings */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-function"
#include <polyscope/polyscope.h>
#pragma GCC diagnostic pop

#include <Eigen/Eigen>
#include <polyscopeHelper.hpp>

Eigen::Vector3f polyscopeScreenCoordsToWorldPos(glm::vec2 screenCoords)
{
    auto view = polyscope::view::getCameraViewMatrix();
    auto proj = polyscope::view::getCameraPerspectiveMatrix();
    glm::vec4 viewport = { 0., 0., polyscope::view::windowWidth, polyscope::view::windowHeight };

    glm::vec3 screenPos3{ screenCoords.x, polyscope::view::windowHeight - screenCoords.y, 0. };
    auto worldPos = glm::unProject(screenPos3, view, proj, viewport);
    return Eigen::Vector3f(worldPos.x, worldPos.y, worldPos.z);
}
