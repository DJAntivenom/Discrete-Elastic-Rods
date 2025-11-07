/**
 * \file DiscreteElasticRod.hpp
 * \brief Class definition of \ref DiscreteElasticRod.
 */

#ifndef __DISCRETE_ELASTIC_ROD_HPP__
#define __DISCRETE_ELASTIC_ROD_HPP__

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
     * \brief Corresponds to n from the paper, i.e. count of all internal vertices.
     */
    uint64_t m_n;

public:
    /**
     * @brief Creates a discrete elastic rod with `n` inner vertices.
     * @param n The number of inner vertices. There will be 2 more vertices,
     *  corresponding to the constrained boundary vertices.
     */
    DiscreteElasticRod(uint64_t n, float theta_zero = 0.f, float theta_n = 0.f);
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
     * @brief Get the stacked vertex velocities.
     * @return A column-vector of size 3*(n+2) with the velocities of the vertices.
    */
    inline const Eigen::Matrix3Xf &getVertexVelocities() const { return m_vertex_velocities; }

    /**
    * @brief Get the stacked edge angles.
    * @return A column-vector of size n with the angles (compared to bishop from) of the vertices.
    */
    inline const Eigen::VectorXf &getEdgeThetas() const { return m_edge_theta; }

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

private:
    void doSymplecticEuler(double delta_time);

    void transportBishopFrame();

    void applyTwist(size_t max_newton_iterations);
};

#endif
