/**
 * \file DiscreteElasticRod.hpp
 * \brief Class definition of \ref DiscreteElasticRod.
 */

#ifndef __POLYSCOPE_HELPER_HPP__
#define __POLYSCOPE_HELPER_HPP__

/* Get rid of annoying unused ... warnings */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-function"
#include <polyscope/polyscope.h>
#pragma GCC diagnostic pop

#include <Eigen/Eigen>

/**
 * \brief Convert screen coordinates to world position
 *
 * Source: polyscope::view::screenCoordsToWorldRay
 * Shortened to return its intermediate value, `worldPos`.
 *
 * Note: This throws warnings on compile. I am not sure why,
 * since this is library code and everything works.
 * This is separate because of that.
 */
Eigen::Vector3f polyscopeScreenCoordsToWorldPos(glm::vec2 screenCoords);

#endif
