#include <DiscreteElasticRod.hpp>

#include <stdexcept>

void DiscreteElasticRod::doSymplecticEuler(double delta_time)
{
    // Calculate forces only at the `n` inner vertices
    Eigen::Matrix3Xf f(3, m_n);
    f.setZero();

    // TODO: Consider if the rod is straight or not (twist branch field m_is_straight_isotropic)
    // TODO: Use the alpha and beta from twist rod

    // Source: [1] 9.1 / Figure 7.
    float alpha = 1.345;
    float beta = 0.789f;

    float start_to_end_theta = m_edge_theta[m_edge_theta.size() - 1] - m_edge_theta[0];

    Eigen::Matrix3Xf kb = getBinormals();
    Eigen::Matrix3Xf tangents = getTangents();

    // m_is_straight_isotropic case:

    for (uint64_t j = 1; j <= m_n; j++)
    {
        // For each node, consider its neighboring vertices
        for (uint64_t i = j - 1; i <= j + 1; i++)
        {
            // Position contribution
            f(j) += -2 * alpha / m_l_i[i];

            // Angle contribution
        }
    }

    (void) delta_time;
    throw std::runtime_error(std::string("DiscreteElasticRod::doSymplecticEuler(double): Not implemented"));
}
