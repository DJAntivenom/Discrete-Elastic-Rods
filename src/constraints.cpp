#include <DiscreteElasticRod.hpp>
#include <Optimization.h>
#include <common.h>

#include <stdexcept>

void DiscreteElasticRod::applyConstraints(size_t max_newton_iterations, double h) {
    using MatType = Eigen::SparseMatrix<Float>;
    const auto C = [this](const VectorX &x_vec)
        {
            const auto x = x_vec.reshaped(3, m_n + 2);
            VectorX c(m_n + 3);
            c[0] = (x.col(0) - m_x0_target).squaredNorm();
            c[m_n + 2] = (x.col(m_n + 1) - m_xnp1_target).squaredNorm();

            c.segment(1, m_n + 1) = (x.rightCols(m_n + 1) - x.leftCols(m_n + 1)).colwise().squaredNorm().transpose()
                - m_edge_length.cwiseSquare();

            return c;
};
    const auto C_partial_x = [this](const VectorX &x_vec)
        {
            const auto x = x_vec.reshaped(3, m_n + 2);

            std::vector<Eigen::Triplet<Float>> triplets;
            {
                Vector3 partialx_0 = 2 * (x.col(0) - m_x0_target);
                Vector3 partialx_np1 = 2 * (x.col(m_n + 1) - m_xnp1_target);
                for (uint8_t i = 0; i < 3u; ++i)
                {
                    triplets.emplace_back(0, i, partialx_0[i]);
                    triplets.emplace_back(m_n + 2, 3 * (m_n + 1) + i, partialx_np1[i]);
                }
            }

            for (uint64_t edge_first_vertex = 0u; edge_first_vertex <= m_n; ++edge_first_vertex)
            {
                Vector3 edge = x.col(edge_first_vertex + 1) - x.col(edge_first_vertex);
                Vector3 partial_first_vertex = -2 * edge;
                Vector3 partial_second_vertex = 2 * edge;

                for (uint8_t i = 0; i < 3; ++i)
                {
                    triplets.emplace_back(1 + edge_first_vertex, 3 * edge_first_vertex + i, partial_first_vertex[i]);
                    triplets.emplace_back(1 + edge_first_vertex, 3 * (edge_first_vertex + 1) + i, partial_second_vertex[i]);
                }
            }

            MatType mat(m_n + 3, 3 * (m_n + 2));
            mat.setFromTriplets(triplets.cbegin(), triplets.cend());
            return mat;
        };
    MatType M(3 * (m_n + 2), 3 * (m_n + 2));
    {
        std::vector<Eigen::Triplet<Float>> triplets;
        for (uint64_t i = 0; i < m_n + 2; ++i)
        {
            triplets.emplace_back(i * 3 + 0, i * 3 + 0, m_vertex_mass[i]);
            triplets.emplace_back(i * 3 + 1, i * 3 + 1, m_vertex_mass[i]);
            triplets.emplace_back(i * 3 + 2, i * 3 + 2, m_vertex_mass[i]);
        }
        M.setFromSortedTriplets(triplets.cbegin(), triplets.cend());
    }

    VectorX x = m_vertex_positions.reshaped(3 * (m_n + 2), 1) + m_vertex_velocities.reshaped(3 * (m_n + 2), 1) * h;
    for (size_t step = 0; step < max_newton_iterations; ++step)
    {
        auto current_residual = C(x);
        if (current_residual.norm() < 1e-8)
        {
            print_debug("Constraints satisfied after " + std::to_string(step) + " steps\n");
            break;
        }

        const auto partial = C_partial_x(x);

        const MatType A(h * h * partial * M.cwiseInverse() * partial.transpose());
        Eigen::SimplicialLDLT<MatType> solver(A);

        VectorX delta_lambda = solver.solve(-current_residual);
        VectorX delta_x = h * h * M.cwiseInverse() * partial.transpose() * delta_lambda;
        x += delta_x;
    }

    m_vertex_velocities = 1 / h * (x.reshaped(3, m_n + 2) - m_vertex_positions);
    m_vertex_positions = x.reshaped(3, m_n + 2);
    assert(!m_vertex_positions.hasNaN());
    assert(m_vertex_positions.allFinite());
}
