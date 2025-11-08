#include <DiscreteElasticRod.hpp>

#include <iostream>

DiscreteElasticRod::DiscreteElasticRod(uint64_t n, float theta_zero, float theta_n) :
    m_vertex_positions(3, (n + 2)),
    m_vertex_velocities(3, (n + 2)),
    m_bishop_frame(3, n + 1),
    m_edge_theta(n + 1),
    m_edge_length(n + 1),
    m_n(n)
{
    Eigen::Vector3f start{ -0.5f, 0.f, 0.f };
    Eigen::Vector3f increment{ 0.1f, 0.0f, 0.0f };
    for (uint64_t i = 0; i < n + 2; ++i)
    {
        m_vertex_positions.col(i) = start + i * increment;
        if (i > 0)
            m_edge_length[i - 1] = 0.1;
    }
    m_bishop_frame_vector = Eigen::Vector3f(0.f, 0.f, 0.1f); //perpendicular to first edge

    m_total_rod_length = m_edge_length.sum();

    m_vertex_velocities.setZero(); //starting at rest

    /**
     * TODO: once we don't just use naturally straigth rods anymore, this has to be done
     * differently.
     */
    m_bishop_frame.colwise() = m_bishop_frame_vector;

    m_edge_theta.setZero();
    m_edge_theta(0) = theta_zero;
    m_edge_theta(n) = theta_n;
}

void DiscreteElasticRod::update(double delta_time, size_t max_newton_iterations)
{
    static float base_time = 0;
    const float y = std::sin(base_time) * 0.1f;
    base_time += delta_time;

    m_vertex_positions.col(m_n / 2).y() = y;

    /* algorithm outline */

    /// 4., 5. apply torque and integrate rigid body
    /// is handled by just setting the handle to some position

    /// 6., 7. compute forces (given in forumla 11 and above (sec. 7.1))
    /// integrate centerline => apply symplectic euler, see ex.1 handout for formula
    /// initial velocity is zero
    doSymplecticEuler(delta_time);

    /// 8. TODO: enforce constraints to guarantee inextensibility

    /// 9. TODO: handle collisions (not discussed in detail in paper, but reference in paper)

    /// 10. update natural bishop frame, i.e. apply rotation (P_i in paper)
    transportBishopFrame();

    /// 11. quasistatic material frame update (Newton according to equation 4 in paper)
    /// => Use newton solver from exercises
    applyTwist(max_newton_iterations);
}

void DiscreteElasticRod::randomizeVertexPositions()
{
    m_vertex_positions.block(1, 0, 2, m_n + 2).setRandom() *= 0.1f;

    transportBishopFrame();
}
