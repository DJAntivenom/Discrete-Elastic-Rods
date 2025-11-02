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
     * \brief Column-vector of stacked vertex positions, size 3*(n+2).
     */
    Eigen::VectorXf m_vertex_positions;

    /*
     * \brief Column-vector of stacked vertex velocities, size 3*(n+2).
    */
    Eigen::VectorXf m_vertex_velocities;

    /*
     * \brief Initial bishop frame vector at edge 0.
    */
    Eigen::Vector3f bishop_frame_vector;

    /*
     * \brief Angle of rotation for each edge compared to bishop frame.
    */
    Eigen::VectorXf m_edge_theta;

    /**
     * \brief Corresponds to (n+2) from the paper, i.e. count of all vertices.
     */
    uint64_t m_num_vertices;

    /*
     * \brief describes if boundary condition is clamped (true) or stressfree (false).    
    */
    bool is_clamped;

public:
    /**
     * @brief Creates a discrete elastic rod with `n` inner vertices.
     * @param n The number of inner vertices. There will be 2 more vertices,
     *  corresponding to the constrained boundary vertices.
     */
    DiscreteElasticRod(uint64_t n, bool clamped = false, float theta_zero = 0.f, float theta_n = 0.f);
    virtual ~DiscreteElasticRod() = default;

    /**
     * \brief Performs one timestep of the simulation.
     * \param delta_time The number of seconds that have passed since the last update.
     */
    void update(double delta_time);

    /**
     * @brief Get the stacked vertex positions.
     * @return A column-vector of size 3*(n+2) with the positions of the vertices.
     */
    inline const Eigen::VectorXf &getVertexPositions() { return m_vertex_positions; }

    /**
     * @brief Get the stacked vertex velocities.
     * @return A column-vector of size 3*(n+2) with the velocities of the vertices.
    */
    inline const Eigen::VectorXf &getVertexVelocities() { return m_vertex_velocities; }

    /**
    * @brief Get the stacked edge angles.
    * @return A column-vector of size n with the angles (compared to bishop from) of the vertices.
    */
    inline const Eigen::VectorXf &getEdgeThetas() { return m_edge_theta; }

private:
    inline auto getVertex(uint64_t vertex_index) { return m_vertex_positions.segment<3>(3 * vertex_index); }
    inline auto getVelocity(uint64_t vertex_index) { return m_vertex_velocities.segment<3>(3 * vertex_index); }
    inline auto getTheta(uint64_t edge_index) { return m_edge_theta(edge_index); }
};

#endif
