#include <DiscreteElasticRod.hpp>
#include <common.h>
#include <stdexcept>

// #define USE_GRAVITY

/**
 * Curvature binormal derivative by other node, $\nabla_i(kb)_i$, from paper chapter 7.1.
 * @return $[\nabla_{i-1}(kb)_i, \nabla_i(kb)_i, \nabla_{i+1}(kb)_i]$
 */
static std::vector<Matrix3> nabla_curvature_binormal(
    const uint64_t i,
    const Matrix3X &edges,
    const Matrix3X &binormals_padded,
    const VectorX &rest_edge_lengths)
{
    const Vector3 e_i = edges.col(i);
    const Vector3 e_i_min_1 = edges.col(i);
    const Vector3 kb_i = binormals_padded.col(i);

    Float denominator = rest_edge_lengths[i-1] * rest_edge_lengths[i] + e_i_min_1.dot(e_i);

    Matrix3 skew_e_i;
    skew_e_i << 0, -e_i[2], e_i[1],
        e_i[2], 0, -e_i[0],
        -e_i[1], e_i[0], 0;
    Matrix3 skew_e_i_min_1;
    skew_e_i_min_1 << 0, -e_i_min_1[2], e_i_min_1[1],
        e_i_min_1[2], 0, -e_i_min_1[0],
        -e_i_min_1[1], e_i_min_1[0], 0;

    Matrix3 nabla_min_1 = (2 * skew_e_i + kb_i * e_i.transpose()) / denominator;
    // TODO: Fix this computation if paper has typo?
    Matrix3 nabla_plus_1 = (2 * skew_e_i_min_1 + kb_i * e_i_min_1.transpose()) / denominator;
    Matrix3 nabla = -(nabla_min_1 + nabla_plus_1);

    return { nabla_min_1, nabla, nabla_plus_1 };
}

void DiscreteElasticRod::doSymplecticEuler(double delta_time)
{
    print_debug("enter");

    // Calculate forces only at the `n` inner vertices
    Matrix3X f(3, m_n + 2);
    f.setZero();

    Float start_to_end_theta = m_edge_theta[m_edge_theta.size() - 1] - m_edge_theta[0];

    Matrix3X binormals = getBinormals();
    Matrix3X edges = getEdges();
    Matrix3X tangents = getTangents();

    VectorX actual_edge_lengths = edges.colwise().norm();
    float actual_rod_length = actual_edge_lengths.sum();

    // Binormals padded with start and end
    Matrix3X binormals_padded(3, m_n + 2);
    binormals_padded.setZero();
    binormals_padded.block(0, 1, 3, m_n) = binormals;

    // TODO: Implement case m_is_straight_isotropic = false

    // --- Calculate Forces
    if (m_is_straight_isotropic || true) // TODO: Undo debug if
    {
        for (uint64_t i = 1; i <= m_n; i++)
        {
            auto nabla_kb = nabla_curvature_binormal(i, edges, binormals_padded, m_edge_length);

            // Equation (9) of paper
            Vector3 nabla_psi_i_min_1 = binormals_padded.col(i) / (2 * m_edge_length[i - 1]);
            Vector3 nabla_psi_i_plus_1 = -binormals_padded.col(i) / (2 * m_edge_length[i]);
            Vector3 nabla_psi_i = -(nabla_psi_i_min_1 + nabla_psi_i_plus_1);
            std::vector<Vector3> nabla_psi = { nabla_psi_i_min_1, nabla_psi_i, nabla_psi_i_plus_1 };

            // For each node, consider its neighboring vertices
            for (uint64_t j = i - 1; j <= i + 1; j++)
            {
                // j to index into nabla vectors with
                const uint64_t j_nabla = j - i + 1;

                // Position contribution
                Vector3 force_position = -2 * m_alpha / m_l_i[i] * nabla_kb[j_nabla].transpose() * binormals_padded.col(i);
                f.col(j) += force_position; // TODO: Is j or i correct?

                // Angle contribution
                Vector3 force_angle = +m_beta * start_to_end_theta / m_total_rod_length * nabla_psi[j_nabla];
                f.col(j) += force_angle; // TODO: Is j or i correct?
            }

            // Gravity, pulling down
#ifdef USE_GRAVITY
            f(1, i) = m_vertex_mass[i] * -9.81f;
#endif

            // TODO: Looks like there's an off by one error, node i==8 has force instead of i==9
            print_debug("Force on node " + std::to_string(i) + ":");
            print_debug(f.col(i));
        }
    }
    else
    {
        // TODO: Implement non-isotropic case
        throw std::runtime_error(std::string(
            "DiscreteElasticRod::doSymplecticEuler(double): Not implemented for non-straight or non-isotropic case"));
    }

    // --- Apply symplectic euler update
    for (uint64_t i = 1; i <= m_n; i++)
    {
        m_vertex_velocities.col(i) += delta_time * f.col(i) / m_vertex_mass[i];
        m_vertex_positions.col(i) += m_vertex_velocities.col(i) * delta_time;
    }

    print_debug("exit");
}
