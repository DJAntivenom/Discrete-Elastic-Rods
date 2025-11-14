/**
 * \file DiscreteElasticRod.hpp
 * \brief Class definition of \ref DiscreteElasticRod.
 */

#ifndef __DISCRETE_ELASTIC_ROD_HPP__
#define __DISCRETE_ELASTIC_ROD_HPP__

#include "Optimization.h"

#include <string>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-function"
#include <polyscope/surface_mesh.h>
#pragma GCC diagnostic pop

#include <Eigen/Eigen>

class DiscreteElasticRod
{
private:
    /**
     * \brief Matrix of vertex positionts, 3x(n+2)
     *
     * Stores positions as column vectors of all vertices (including boundaries).
     */
    Eigen::Matrix3Xf m_vertex_positions;

    /**
     * \brief Matrix of vertex velocities, 3x(n+2)
     *
     * Stores velocities as column vectors of all vertices (including boundaries).
     */
    Eigen::Matrix3Xf m_vertex_velocities;

    /**
     * \brief Initial bishop frame vector u at edge 0.
     */
    Eigen::Vector3f m_bishop_frame_vector;

    /**
     * \brief Column i is bishop-frame vector u_i.
     *
     * Dimension: 3 x (n+1)
     *
     * The frames are updated after every update step in `transportBishopFrame`.
     */
    Eigen::Matrix3Xf m_bishop_frame;

    /**
     * \brief Angle of rotation theta for each edge compared to bishop frame.
     *
     * Dimension: (n+1)
     */
    Eigen::VectorXf m_edge_theta;

    /**
     * \brief Length of each edge, calculated in constructor
     *
     * Dimension: (n+1)
     */
    Eigen::VectorXf m_edge_length;

    /**
     * \brief Half the size of the voronoi region associated with each vertex.
     *
     * Dimension: (n + 2) but not useful at boundaries
     */
    Eigen::VectorXf m_l_i;

    /**
     * \brief The total length of the rod.
     */
    float m_total_rod_length;

    /**
     * \brief Corresponds to n from the paper, i.e. count of all internal vertices.
     */
    uint64_t m_n;

    /**
     * \brief Radius of the rod.
     */
    float m_radius;

    /**
     * @brief Set when positions of vertices are randomized or moved.
     */
    bool m_is_straight_isotropic;

    /**
     * \brief Material parameters.
     */
    float m_alpha, m_beta;

    /**
     * \brief \overbar{B}^j from the paper.
     */
    Eigen::Matrix2f m_B_matrix;

    /**
     * @brief \overbar{\omega} from the paper.
     */
    Eigen::Matrix4Xf m_w_overbar;

public:
    /**
     * @brief Creates a discrete elastic rod with `n` inner vertices.
     * @param n The number of inner vertices. There will be 2 more vertices,
     *  corresponding to the constrained boundary vertices.
     */
    DiscreteElasticRod(uint64_t n,
                       float alpha, float beta,
                       float radius = 0.5f, float theta_zero = 0.f, float theta_n = 0.f);
    virtual ~DiscreteElasticRod() = default;

    /**
     * \brief Performs one timestep of the simulation.
     * \param delta_time The number of seconds that have passed since the last update.
     */
    void update(double delta_time, size_t max_newton_iterations);

    /**
     * @brief Get the stacked vertex positions.
     * @return A column-vector of size 3*(n+2) with the positions of the vertices.
     */
    inline const Eigen::Matrix3Xf &getVertexPositions() const { return m_vertex_positions; }

    /**
     * @brief Set vertex position.
     */
    void setVertexPosition(uint64_t vertex_index, const Eigen::Vector3f &new_position);

    /**
     * @brief Get the vertex velocities.
     * @return A matrix of size 3x(n+2) with the velocities of the vertices.
     */
    inline const Eigen::Matrix3Xf &getVertexVelocities() const { return m_vertex_velocities; }

    /**
     * \brief Get u_i in the columns of the returned matrix.
     */
    inline const Eigen::Matrix3Xf &getBishopFrame() const { return m_bishop_frame; }

    /**
    * @brief Get the stacked edge angles.
    * @return A column-vector of size n with the angles (compared to bishop from) of the vertices.
    */
    inline const Eigen::VectorXf &getEdgeThetas() const { return m_edge_theta; }

    /**
    * @brief Set the edge angle boundary conditions.
    */
    inline void setBoundaryEdgeThetas(float theta_start, float theta_end)
    {
        m_edge_theta[0] = theta_start;
        m_edge_theta[m_edge_theta.size() - 1] = theta_end;
    }

    /**
     * \brief Get the lengths of the edges.
     */
    inline const Eigen::VectorXf &getEdgeLengths() const { return m_edge_length; }

    /**
     * @brief Get $e_i$ from the paper, i.e. the vectors representing segments between vertices.
     * @return A matrix, where each of the (n+1) columns is one edge.
     */
    inline Eigen::Matrix3Xf getEdges() const
    {
        return m_vertex_positions.block(0, 1, 3, m_n + 1)
            - m_vertex_positions.block(0, 0, 3, m_n + 1);
    }

    /**
     * \brief Get $t_i$ from the paper, i.e. discrete tangents.
     * \return A matrix, where each of the (n+1) columns is one tangent.
     */
    inline Eigen::Matrix3Xf getTangents() const
    {
        return getEdges().colwise().normalized();
    }

    /**
     * @brief Get kappa_i from the paper, i.e. discrete curvature at vertices.
     * @return A vector of length n containing discrete curvature at each inner vertex.
     */
    Eigen::VectorXf getCurvature() const
    {
        const auto edges = getTangents();

        /* there is no per-column dot product function => do dot-product "by hand" */
        const auto angles = (edges.block(0, 0, 3, m_n).array() *
                             edges.block(0, 1, 3, m_n).array()).colwise().sum().acos();

        /* don't forget to go back to matrices, like I initially did ^^ */
        return 2. * (angles * 0.5).tan().matrix();
    }

    /**
     * @brief Get (kb)_i from the paper, i.e. discrete curvature binormal.
     * @return A 3 x n matrix, where column i is the binormal at inner vertex i, i.e.
     * global vertex (i+1), because our boundary vertex is 0.
     */
    Eigen::Matrix3Xf getBinormals() const {
        /* instead of normalizing the corss-product and scaling by curvature,
           we use (1) from the paper which is robust if edges are collinear
           (normalizing of a 0-vector, i.e. result of collinear x-product, is
           undefined)*/
        const auto edges = getEdges();
        Eigen::Matrix3Xf binormals(3, m_n);
        for (uint32_t col_index = 0; col_index < m_n; ++col_index)
        {
            binormals.col(col_index) = 2 * edges.col(col_index).cross(edges.col(col_index + 1));
            binormals.col(col_index) *= 1 / (edges.col(col_index).norm() * edges.col(col_index + 1).norm() +
                                             edges.col(col_index).dot(edges.col(col_index + 1)));
        }
        return binormals;
    }

    /**
     * \brief Get omega_i from the paper
     * \return A 4 x (n + 2) matrix with omega at the vertices, where every column is
     * {omega_i^[i-1].x, omega_i^[i-1].y, omega_i^[i].x, omega_i^[i].y}.
     *
     * The values at the boundary vertices is undefined.
     */
    Eigen::Matrix4Xf getMaterialCurvature() const;

    /**
     * \brief Get the bishop frame {x, u, v} at alpha along the rod.
     *
     * The frame is represented in reduced form as a 3x3 matrix, where
     * the first row is x, the second is u and the third row is v.
     * \param alpha The ratio describing how far along the rod the frame
     * should be interpolated.
     * \return The position x interpolated along the rod in the first row
     * followed by u and v.
     */
    Eigen::Matrix3f getInterpolatedBishopFrame(double alpha);

    /**
     * \brief Randomize the positions of the vertices within an area around the initial position.
     */
    void randomizeVertexPositions();

    /**
     * @brief Whether this rod is straight and isotropic when at rest.
     */
    inline bool is_straight_isotropic() const { return m_is_straight_isotropic; }

    /**
     * @brief Create a surface mesh from the current positions and edge thetas
     * @param name This is passed to polyscope for later retrieval.
     * @param vertices_per_ring How smooth the rod should be drawn.
     */
    polyscope::SurfaceMesh *registerSurfaceMesh(const std::string &name,
                                                uint32_t vertices_per_ring = 8u) const;

private:
    void doSymplecticEuler(double delta_time);

    void transportBishopFrame();

    void applyTwist(size_t max_newton_iterations);

    bool twistEnergy(const Optimization::VectorXf &theta, double &energy) const;
    bool twistGradient(const Optimization::VectorXf &theta, Optimization::VectorXf &gradient) const;
    bool twistHessian(const Optimization::VectorXf &theta, Optimization::TripletListF &hessian) const;
};

#endif
