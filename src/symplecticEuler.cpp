#include <DiscreteElasticRod.hpp>
#include <common.h>
#include <stdexcept>

/**
 * Curvature binormal derivative by other node, $\nabla_i(kb)_i$, from paper chapter 7.1.
 * @return $[\nabla_{i-1}(kb)_i, \nabla_i(kb)_i, \nabla_{i+1}(kb)_i]$
 */
static std::vector<Eigen::Matrix3f> nabla_curvature_binormal(
    const uint64_t i,
    const Eigen::Matrix3Xf &edges,
    const Eigen::Matrix3Xf &binormals_padded)
{
    const Eigen::Vector3f e_i = edges.col(i);
    const Eigen::Vector3f e_i_min_1 = edges.col(i);
    const Eigen::Vector3f kb_i = binormals_padded.col(i);

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
    print_debug("[doSymplecticEuler] enter");

    using namespace Eigen;

    // Calculate forces only at the `n` inner vertices
    Matrix3Xf f(3, m_n+2);
    f.setZero();

    float start_to_end_theta = m_edge_theta[m_edge_theta.size() - 1] - m_edge_theta[0];

    print_debug("[doSymplecticEuler] get rod values");
    Matrix3Xf binormals = getBinormals();
    Matrix3Xf edges = getEdges();
    Matrix3Xf tangents = getTangents();

    // Binormals padded with start and end
    Matrix3Xf binormals_padded(3, m_n + 2);
    binormals_padded.setZero();
    binormals_padded.block(0, 1, 3, m_n) = binormals;

    // TODO: Implement case m_is_straight_isotropic = false

    // --- Calculate Forces
    print_debug("[doSymplecticEuler] calculate forces");

    if (m_is_straight_isotropic || true) // TODO: Undo debug if
    {
        for (uint64_t i = 1; i <= m_n; i++)
        {
            print_debug("[doSymplecticEuler] nabla kb");
            auto nabla_kb = nabla_curvature_binormal(i, edges, binormals_padded);

            // Equation (9) of paper
            print_debug("[doSymplecticEuler] nabla psi");
            Vector3f nabla_psi_i_min_1 = binormals_padded.col(i) / (2 * m_edge_length[i - 1]);
            Vector3f nabla_psi_i_plus_1 = binormals_padded.col(i) / (2 * m_edge_length[i]);
            Vector3f nabla_psi_i = -(nabla_psi_i_min_1 + nabla_psi_i_plus_1);
            std::vector<Vector3f> nabla_psi = {nabla_psi_i_min_1, nabla_psi_i, nabla_psi_i_plus_1};

            print_debug("[doSymplecticEuler] get node force: " + std::to_string(i));
            // For each node, consider its neighboring vertices
            for (uint64_t j = i - 1; j <= i + 1; j++)
            {
                // j to index into nabla vectors with
                const uint64_t j_nabla = j - i + 1;

                // Position contribution
                print_debug("[doSymplecticEuler] get node force - 1");
                (void)m_l_i[j];
                print_debug("[doSymplecticEuler] get node force - 2");
                (void)nabla_kb[j_nabla].transpose();
                print_debug("[doSymplecticEuler] get node force - 3");
                f.col(i) = Vector3f::Zero();
                print_debug("[doSymplecticEuler] get node force - position");
                f.col(i) += -2 * m_alpha / m_l_i[j] * nabla_kb[j_nabla].transpose() * binormals_padded.col(j);

                // Angle contribution
                print_debug("[doSymplecticEuler] get node force - angle");
                f.col(i) += m_beta * start_to_end_theta / m_total_rod_length * nabla_psi[j_nabla];
            }
        }
    }
    else
    {
        // TODO: Implement non-isotropic case
        throw std::runtime_error(std::string(
            "DiscreteElasticRod::doSymplecticEuler(double): Not implemented for non-straight or non-isotropic case"));
    }

    // --- Apply symplectic euler update
    print_debug("[doSymplecticEuler] apply euler update");
    for (uint64_t i = 1; i <= m_n; i++)
    {
        print_debug("[doSymplecticEuler] apply euler update: " + std::to_string(i));

        print_debug("[doSymplecticEuler] apply euler update - velocity");
        m_vertex_velocities.col(i) += delta_time * f.col(i) / m_vertex_mass[i];
        print_debug("[doSymplecticEuler] apply euler update - position");
        m_vertex_positions.col(i) += m_vertex_velocities.col(i) * delta_time;
    }

    print_debug("[doSymplecticEuler] exit");
}
