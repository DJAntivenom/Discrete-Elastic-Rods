#ifndef DISCRETE_ELASTIC_RODS_COMMON_H
#define DISCRETE_ELASTIC_RODS_COMMON_H

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-function"
#include <polyscope/surface_mesh.h>
#pragma GCC diagnostic pop

#include <Eigen/Eigen>

#ifndef NDEBUG
#define PRINT_DEBUG
#endif

#ifdef PRINT_DEBUG
#define print_debug(str) std::cout << "[" << __func__ << "] " << (str) << std::endl
#else
#define print_debug(str) \
    do                   \
    {                    \
    } while (0)
#endif

/////////////
/// TYPES ///
/////////////

/**
 * @brief The datatype of the values stored in this rod.
 */
using Float = double;

/// Vectors
using Vector2 = Eigen::Vector2<Float>;
using Vector3 = Eigen::Vector3<Float>;
using VectorX = Eigen::VectorX<Float>;

/// Fixed-size matrices
using Matrix2 = Eigen::Matrix2<Float>;
using Matrix3 = Eigen::Matrix3<Float>;

/// Partially-dynamic matrices
using Matrix3X = Eigen::Matrix3X<Float>;
using MatrixX3 = Eigen::MatrixX3<Float>;
using Matrix4X = Eigen::Matrix4X<Float>;

/// Rotations
using AngleAxis = Eigen::AngleAxis<Float>;
using AffineTransform3D = Eigen::Transform<Float, 3, Eigen::Affine>;
using GeneralTransform3D = Eigen::Transform<Float, 3, Eigen::Projective>;
using IsometricTransform3D = Eigen::Transform<Float, 3, Eigen::Isometry>;

#endif // DISCRETE_ELASTIC_RODS_COMMON_H
