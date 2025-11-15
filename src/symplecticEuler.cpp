#include <DiscreteElasticRod.hpp>

#include <stdexcept>

/**
 * Curvature binormal derivative by other node, $\nabla_i(kb)_i$, from paper chapter 7.1.
 * @return $[\nabla_{i-1}(kb)_i, \nabla_i(kb)_i, \nabla_{i+1}(kb)_i]$
 */
static std::vector<Eigen::Matrix3f> nabla_curvature_binormal(
    const uint64_t i,
    const Eigen::Matrix3Xf &edges,
    const Eigen::Matrix3Xf &binormals)
{
    const Eigen::Vector3f e_i = edges.col(i);
    const Eigen::Vector3f e_i_min_1 = edges.col(i);
    const Eigen::Vector3f kb_i = binormals.col(i);

    float denominator = e_i_min_1.norm() * e_i.norm() + e_i_min_1.dot(e_i);

    Eigen::Matrix3f skew_e_i;
    skew_e_i << 0, -e_i[2], e_i[1],
            e_i[2], 0, -e_i[0],
            -e_i[1], e_i[0], 0;
    Eigen::Matrix3f skew_e_i_min_1;
    skew_e_i_min_1 << 0, -e_i_min_1[2], e_i_min_1[1],
            e_i_min_1[2], 0, -e_i_min_1[0],
            -e_i_min_1[1], e_i_min_1[0], 0;

    Eigen::Matrix3f nabla_min_1 = 2 * skew_e_i + kb_i * e_i.transpose() / denominator;
    // TODO: Fix this computation if paper has typo?
    Eigen::Matrix3f nabla_plus_1 = 2 * skew_e_i_min_1 + kb_i * e_i_min_1.transpose() / denominator;
    Eigen::Matrix3f nabla = -(nabla_min_1 + nabla_plus_1);

    return {nabla_min_1, nabla, nabla_plus_1};
}

void DiscreteElasticRod::doSymplecticEuler(double delta_time)
{
    using namespace Eigen;

    // Calculate forces only at the `n` inner vertices
    Matrix3Xf f(3, m_n);
    f.setZero();

    // TODO: Consider if the rod is straight or not (twist branch field m_is_straight_isotropic)

    const float start_to_end_theta = m_edge_theta[m_edge_theta.size() - 1] - m_edge_theta[0];

    const Matrix3Xf kb = getBinormals();
    const Matrix3f edges = getEdges();
    const Matrix3Xf tangents = getTangents();

    const float L = m_edge_length.sum();

    // TODO: Implement case m_is_straight_isotropic = false

    // --- Calculate Forces

    if (m_is_straight_isotropic)
    {
        for (uint64_t i = 1; i <= m_n; i++)
        {
            auto nabla_kb = nabla_curvature_binormal(i, edges, kb);

            // Equation (9) of paper
            const Vector3f nabla_psi_i_min_1 = kb.col(i) / (2 * m_edge_length[i - 1]);
            const Vector3f nabla_psi_i_plus_1 = kb.col(i) / (2 * m_edge_length[i]);
            const Vector3f nabla_psi_i = -(nabla_psi_i_min_1 + nabla_psi_i_plus_1);
            std::vector<Vector3f> nabla_psi = {nabla_psi_i_min_1, nabla_psi_i, nabla_psi_i_plus_1};

            // For each node, consider its neighboring vertices
            for (uint64_t j = i - 1; j <= i + 1; j++)
            {
                // j to index into nabla vectors with
                const uint64_t j_nabla = j - i + 1;

                // Position contribution
                f(i, placeholders::all) += -2 * m_alpha / m_l_i[j] * nabla_kb[j_nabla].transpose() * kb.col(j);

                // Angle contribution
                f(i, placeholders::all) += m_beta * start_to_end_theta / L * nabla_psi[j_nabla];
            }
        }
    }
    else
    {
        throw std::runtime_error(std::string(
            "DiscreteElasticRod::doSymplecticEuler(double): Not implemented for non-straight or non-isotropic case"));
    }

    // --- Apply update


    (void) delta_time;
    throw std::runtime_error(std::string("DiscreteElasticRod::doSymplecticEuler(double): Not implemented"));
}
