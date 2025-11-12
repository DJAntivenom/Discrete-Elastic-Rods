#include <DiscreteElasticRod.hpp>

#include <stdexcept>

void DiscreteElasticRod::applyTwist(size_t max_newton_iterations)
{
    if (m_is_straight_isotropic)
    {
        /** TODO: Once we have different boundary conditions we need to calculate the twist differently for those */

        /** theta_i = c * l_i + theta^(i - 1), where c = (theta_n - theta_0) / 2L */
        const float c = (m_edge_theta[m_n] - m_edge_theta[0]) / (2 * m_total_rod_length);
        for (uint32_t i = 1; i < m_n; ++i)
        {
            m_edge_theta[i] = c * m_l_i[i] + m_edge_theta[i - 1];
        }
    }
    else
    {
        (void)max_newton_iterations;
        throw std::runtime_error(std::string("DiscreteElasticRod::applyTwist(): Not implemented for non-straight rods"));
    }
}
