/**
 * \file DiscreteElasticRod.hpp
 * \brief Class definition of \ref DiscreteElasticRod.
 */

#ifndef __DISCRETE_ELASTIC_ROD_HPP__
#define __DISCRETE_ELASTIC_ROD_HPP__

#include "common.h"
#include "Optimization.h"
#include "InitialConfiguration.hpp"

#include <string>

class DiscreteElasticRod
{
public:
    /// Optimizer
    using Opt = Optimization<Float>;

private:
    DiscreteElasticRod(Matrix3X &vertex_positions,
                                        VectorX &edge_lengths,
                                        Matrix3X &vertex_velocities,
                                        Vector3 &bishop_frame_vector,
                                        Matrix3X &bishop_frame,
                                        VectorX &edge_theta,
                                        VectorX &vertex_mass,
                                        VectorX &l_i,
                                        Matrix2 &B_matrix,
                                        Matrix4X &w_overbar,
                                        int n,
                                        Float radius,
                                        bool is_straight_isotropic,
                                        Float alpha,
                                        Float beta);
    /**
     * \brief Matrix of vertex positionts, 3x(n+2)
     *
     * Stores positions as column vectors of all vertices (including boundaries).
     */
    Matrix3X m_vertex_positions;

    /**
     * \brief Matrix of vertex velocities, 3x(n+2)
     *
     * Stores velocities as column vectors of all vertices (including boundaries).
     */
    Matrix3X m_vertex_velocities;

    /**
     * \brief Initial bishop frame vector u at edge 0.
     */
    Vector3 m_bishop_frame_vector;

    /**
     * \brief Column i is bishop-frame vector u_i.
     *
     * Dimension: 3 x (n+1)
     *
     * The frames are updated after every update step in `transportBishopFrame`.
     */
    Matrix3X m_bishop_frame;

    /**
     * \brief Angle of rotation theta for each edge compared to bishop frame.
     *
     * Dimension: (n+1)
     */
    VectorX m_edge_theta;

    /**
     * \brief Length of each edge, calculated in constructor
     *
     * Dimension: (n+1)
     */
    VectorX m_edge_length;

    /**
     * \brief Mass of each vertex
     *
     * Dimension: (n+2)
     */
    VectorX m_vertex_mass;

    /**
     * \brief Half the size of the voronoi region associated with each vertex.
     *
     * Dimension: (n + 2) but not useful at boundaries
     */
    VectorX m_l_i;

    /**
     * \brief The total length of the rod.
     */
    Float m_total_rod_length;

    /**
     * \brief Corresponds to n from the paper, i.e. count of all internal vertices.
     */
    const uint64_t m_n;

    /**
     * \brief Radius of the rod.
     */
    const float m_radius;

    /**
     * @brief Set when positions of vertices are randomized or moved.
     */
    bool m_is_straight_isotropic;

    /**
     * \brief Material parameters.
     */
    const Float m_alpha, m_beta;

    /**
     * \brief \overbar{B}^j from the paper.
     */
    Matrix2 m_B_matrix;

    /**
     * @brief \overbar{\omega} from the paper.
     */
    Matrix4X m_w_overbar;

public:
    /**
     * @brief Creates a discrete elastic rod with `n` inner vertices.
     * @param n The number of inner vertices. There will be 2 more vertices,
     *  corresponding to the constrained boundary vertices.
     */
    DiscreteElasticRod(const InitialConfiguration &initial_configuration,
                       Float alpha, Float beta, Float total_mass = 0.1);

    virtual ~DiscreteElasticRod() = default;

    /**
     * \brief Performs one timestep of the simulation.
     * \param delta_time The number of seconds that have passed since the last update.
     */
    void update(double delta_time, size_t max_newton_iterations);

    /*
     * \brief the following two methods split the simulation time_step in half, before and after the time_step respectively.
     */

    void preCollisionUpdate(double delta_time, size_t max_newton_iterations);

    void calculateContactForces(const int max_iter, DiscreteElasticRod &other, const bool same, const Float delta_time);

    void postCollisionUpdate(double delta_time, size_t max_newton_iterations);

    inline uint64_t getN() const { return m_n; }

    /**
     * @brief Get the stacked vertex positions.
     * @return A column-vector of size 3*(n+2) with the positions of the vertices.
     */
    inline const Matrix3X &getVertexPositions() const { return m_vertex_positions; }

    /**
     * @brief Set vertex position.
     */
    void setVertexPosition(uint64_t vertex_index, const Vector3 &new_position);

    /**
     * @brief Get the vertex velocities.
     * @return A matrix of size 3x(n+2) with the velocities of the vertices.
     */
    inline const Matrix3X &getVertexVelocities() const { return m_vertex_velocities; }

    /**
     * \brief Get u_i in the columns of the returned matrix.
     */
    inline const Matrix3X &getBishopFrame() const { return m_bishop_frame; }

    /**
    * @brief Get the stacked edge angles.
    * @return A column-vector of size n with the angles (compared to bishop from) of the vertices.
    */
    inline const VectorX &getEdgeThetas() const { return m_edge_theta; }

    /**
    * @brief Set the edge angle boundary conditions.
    */
    inline void setBoundaryEdgeThetas(Float theta_start, Float theta_end)
    {
        m_edge_theta[0] = theta_start;
        m_edge_theta[m_edge_theta.size() - 1] = theta_end;
    }

    /**
     * \brief Get the lengths of the edges.
     */
    inline const VectorX &getEdgeLengths() const { return m_edge_length; }

    /**
     * @brief Get $e_i$ from the paper, i.e. the vectors representing segments between vertices.
     * @return A matrix, where each of the (n+1) columns is one edge.
     */
    inline Matrix3X getEdges() const
    {
        return m_vertex_positions.block(0, 1, 3, m_n + 1)
            - m_vertex_positions.block(0, 0, 3, m_n + 1);
    }

    /**
     * \brief Get $t_i$ from the paper, i.e. discrete tangents.
     * \return A matrix, where each of the (n+1) columns is one tangent.
     */
    inline Matrix3X getTangents() const
    {
        return getEdges().colwise().normalized();
    }

    inline Float getRadius() const { return m_radius; }

    /**
     * @brief Get kappa_i from the paper, i.e. discrete curvature at vertices.
     * @return A vector of length n containing discrete curvature at each inner vertex.
     */
    VectorX getCurvature() const
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
    Matrix3X getBinormals() const {
        /* instead of normalizing the corss-product and scaling by curvature,
           we use (1) from the paper which is robust if edges are collinear
           (normalizing of a 0-vector, i.e. result of collinear x-product, is
           undefined)*/
        const auto edges = getEdges();
        Matrix3X binormals(3, m_n);
        for (uint32_t col_index = 0; col_index < m_n; ++col_index)
        {
            binormals.col(col_index) = edges.col(col_index).cross(2 * edges.col(col_index + 1));
            binormals.col(col_index) *= 1 / (m_edge_length(col_index) * m_edge_length(col_index + 1) +
                                             edges.col(col_index).dot(edges.col(col_index + 1)));
        }
        return binormals;
    }

    /**
     * \brief Get omega_i from the paper given edge-wise rotations theta
     * \param theta A vector of edge-wise rotations with same dimension as m_edge_theta.
     * \return A 4 x (n + 2) matrix with omega at the vertices, where every column is
     * {omega_i^[i-1].x, omega_i^[i-1].y, omega_i^[i].x, omega_i^[i].y}.
     *
     * The values at the boundary vertices is undefined.
     */
    Matrix4X getMaterialCurvature(const VectorX &theta) const;

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
    Matrix3 getInterpolatedBishopFrame(double alpha);

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

    /**
     * @brief Cuts the current rod into two separate rods. It is expected that you delete the current rod afterwards.
     * @param i The vertex at which the rod should be cut. This cannot be 0 or m_n + 1
     * @return a pair of DiscreteElasticRods
     */
    std::pair<DiscreteElasticRod, DiscreteElasticRod> cutAtVertex(int i);


private:
    void doSymplecticEuler(double delta_time);

    void transportBishopFrame();

    void applyTwist(size_t max_newton_iterations);

    void applyConstraints(size_t max_newton_iterations);

    bool getConstraints(const Opt::VectorX &q_r_x, Opt::Float &energy) const;
    bool getConstraintGradient(const Opt::VectorX &q_r_x, Opt::VectorX &gradient) const;
    bool getConstraintHessian(const Opt::VectorX &q_r_x, Opt::TripletList &hessian) const;

};

#endif
