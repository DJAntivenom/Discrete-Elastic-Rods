#include <DiscreteElasticRod.hpp>
#include <Optimization.h>

#include <cassert>
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
        auto twistEnergy = [this](const Opt::VectorX &theta, Opt::Float &energy)
            {
                /**
                 * TODO: if we have a stressfree boundary this needs to be theta.size() == m_n
                 * If we have two stressfree ends this is not defined? Because then there shouldn't
                 * be any twist, I guess.
                 */
                assert(theta.size() == m_edge_theta.size());


                const auto inv_l_i = m_l_i.segment(1, theta.size() - 1).cwiseInverse();
                const auto m_j_square_div_l =
                    (theta.tail(theta.size() - 1) -
                     theta.head(theta.size() - 1))
                    .cwiseSquare()
                    .cwiseProduct(inv_l_i);

                energy = m_beta * m_j_square_div_l.sum();

                return true;
            };

        auto twistGradient = [this](const Opt::VectorX &theta, Opt::VectorX &gradient)
            {
                /**
                 * TODO: if we have a stressfree boundary this needs to be theta.size() == m_n
                 * If we have two stressfree ends this is not defined? Because then there shouldn't
                 * be any twist, I guess.
                 */
                assert(theta.size() == m_edge_theta.size());

                gradient.setZero(theta.size());

                const VectorX inv_l_i = m_l_i.cwiseInverse();
                const Matrix4X omega = getMaterialCurvature(theta);
                const auto rotation = Eigen::Rotation2D<Float>(M_PI_2);
                const auto rotated_B = rotation * m_B_matrix;

                const VectorX m_j_div_l =
                    (theta.tail(theta.size() - 1) -
                     theta.head(theta.size() - 1))
                    .cwiseProduct(inv_l_i.segment(1, theta.size() - 1));

                /** TODO: adapt to stressfree boundary condition */
                for (uint32_t j = 1; j < m_n; ++j)
                {
                    /** partial W_j / partial theta^j */
                    {
                        const auto omega_j_j = omega.block<2, 1>(2, j);
                        const auto omega_bar_j_j = m_w_overbar.block<2, 1>(2, j);
                        gradient[j] = inv_l_i[j] * omega_j_j.transpose() * rotated_B * (omega_j_j - omega_bar_j_j);
                    }

                    /** partial W_{j+1} / partial theta^j */
                    {
                        const auto omega_jp1_j = omega.block<2, 1>(0, j + 1);
                        const auto omega_bar_jp1_j = m_w_overbar.block<2, 1>(0, j + 1);
                        gradient[j] += inv_l_i[j + 1] * omega_jp1_j.transpose() * rotated_B * (omega_jp1_j - omega_bar_jp1_j);
                    }

                    /** NOTE: need to subtract 1 because m_j[0] = theta^1 - theta^0 */
                    gradient[j] += 2 * m_beta * (m_j_div_l[j - 1] - m_j_div_l[j]);

                    assert(!std::isnan(gradient[j]));
                    assert(std::isfinite(gradient[j]));
                }

                return true;
            };

        auto twistHessian = [this](const Opt::VectorX &theta, Opt::TripletList &hessian)
            {
                /**
                 * TODO: if we have a stressfree boundary this needs to be theta.size() == m_n
                 * If we have two stressfree ends this is not defined? Because then there shouldn't
                 * be any twist, I guess.
                 */
                assert(theta.size() == m_edge_theta.size());

                const VectorX inv_l_i = m_l_i.cwiseInverse();
                const Matrix4X omega = getMaterialCurvature(theta);
                const auto rotation = Eigen::Rotation2D<Float>(M_PI_2);
                const auto rotated_B = rotation.matrix().transpose() * m_B_matrix * rotation;
                assert(!inv_l_i.hasNaN());
                assert(inv_l_i.allFinite());

                /** hessian.emplace_back(row, col, value) */
                /** hessian has dimensions theta.rows() x theta.rows() */
                for (uint32_t j = 1; j < m_n; ++j)
                {
                    hessian.emplace_back(j, j - 1, -2.f * m_beta * inv_l_i[j]);
                    assert(!std::isnan(hessian.back().value()));
                    assert(std::isfinite(hessian.back().value()));
                    hessian.emplace_back(j, j + 1, -2.f * m_beta * inv_l_i[j + 1]);
                    assert(!std::isnan(hessian.back().value()));
                    assert(std::isfinite(hessian.back().value()));

                    float value = 2.f * m_beta * (inv_l_i[j] + inv_l_i[j + 1]);
                    /** partial^2 W_j / (partial theta^j)^2 */
                    {
                        const auto omega_j_j = omega.block<2, 1>(2, j);
                        const auto omega_bar_j_j = m_w_overbar.block<2, 1>(2, j);
                        value += inv_l_i[j] * omega_j_j.transpose() * rotated_B * omega_j_j;
                        value -= inv_l_i[j] * omega_j_j.transpose() * m_B_matrix * (omega_j_j - omega_bar_j_j);
                    }

                    /** partial^2 W_{j+1} / (partial theta^j)^2 */
                    {
                        const auto omega_jp1_j = omega.block<2, 1>(0, j + 1);
                        const auto omega_bar_jp1_j = m_w_overbar.block<2, 1>(0, j + 1);
                        value += inv_l_i[j + 1] * omega_jp1_j.transpose() * rotated_B * omega_jp1_j;
                        value -= inv_l_i[j + 1] * omega_jp1_j.transpose() * m_B_matrix * (omega_jp1_j - omega_bar_jp1_j);
                    }
                    assert(!std::isnan(value));
                    assert(std::isfinite(value));
                    hessian.emplace_back(j, j, value);
                }

                return true;
            };

        Opt opt;
        opt.optimizer = Opt::Optimizer::NEWTON;
        opt.tolerance_exponent = -10;
        opt.objective_function = [&](const Opt::VectorX &theta, Opt::Float &energy)
            {
                return twistEnergy(theta, energy);
            };
        opt.gradient_function = [&](const Opt::VectorX &theta,
                                    Opt::Float &energy,
                                    Opt::VectorX &gradient)
            {
                twistEnergy(theta, energy);
                return twistGradient(theta, gradient);
            };
        opt.hessian_function = [&](const Opt::VectorX &theta,
                                   Opt::Float &energy,
                                   Opt::VectorX &gradient,
                                   Opt::TripletList &hessian)
            {
                twistEnergy(theta, energy);
                twistGradient(theta, gradient);
                return twistHessian(theta, hessian);
            };

        for (size_t step = 0; step < max_newton_iterations; ++step)
        {
            /** TODO: If we have different boundary conditions, we need to only pass the varying values here */
            VectorX theta = m_edge_theta;
            auto result = opt.step(theta);
            switch (result)
            {
            case Opt::OptimizationStatus::SUCCESS:
                m_edge_theta = theta;
                break;
            case Opt::OptimizationStatus::FAILURE:
                std::cout << "Material frame calculation failed" << std::endl;
                return;
            case Opt::OptimizationStatus::CONVERGED:
                return;
            }
        }

        std::cout << "Reached max iterations in material frame optimization" << std::endl;
    }
}
