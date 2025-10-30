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

    /**
     * \brief Corresponds to (n+2) from the paper, i.e. count of all vertices.
     */
    uint64_t m_num_vertices;

public:
    /**
     * @brief Creates a discrete elastic rod with `n` inner vertices.
     * @param n The number of inner vertices. There will be 2 more vertices,
     *  corresponding to the constrained boundary vertices.
     */
    DiscreteElasticRod(uint64_t n);
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

private:
    inline auto getVertex(uint64_t vertex_index) { return m_vertex_positions.segment<3>(3 * vertex_index); }
};

#endif
