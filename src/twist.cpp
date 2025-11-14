#include <DiscreteElasticRod.hpp>
#include <Optimization.h>

#include <cassert>
#include <stdexcept>

bool DiscreteElasticRod::twistEnergy(const Optimization::VectorXf &theta, double &energy)const
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
}

bool DiscreteElasticRod::twistGradient(const Optimization::VectorXf &theta, Optimization::VectorXf &gradient) const
{
    /**
     * TODO: if we have a stressfree boundary this needs to be theta.size() == m_n
     * If we have two stressfree ends this is not defined? Because then there shouldn't
     * be any twist, I guess.
     */
    assert(theta.size() == m_edge_theta.size());

    gradient.setZero(theta.size());

    const Eigen::VectorXf inv_l_i = m_l_i.cwiseInverse();
    const Eigen::Matrix4Xf omega = getMaterialCurvature(theta);
    const auto rotation = Eigen::Rotation2Df(M_PI_2);
    const auto rotated_B = rotation * m_B_matrix;

    const Eigen::VectorXf m_j_div_l =
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
    }

    return true;
}

bool DiscreteElasticRod::twistHessian(const Optimization::VectorXf &theta, Optimization::TripletListF &hessian) const
{
    (void)theta;
    (void)hessian;
    throw std::runtime_error(std::string(__func__) + ": Not implemented");
}

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
        Optimization opt;
        opt.optimizer = Optimization::Optimizer::NEWTON;
        opt.objective_function = [&](const Eigen::VectorXf &theta, double &energy)
            {
                return twistEnergy(theta, energy);
            };
        opt.gradient_function = [&](const Eigen::VectorXf &theta,
                                    double &energy,
                                    Eigen::VectorXf &gradient)
            {
                twistEnergy(theta, energy);
                return twistGradient(theta, gradient);
            };
        opt.hessian_function = [&](const Optimization::VectorXf &theta,
                                   double &energy,
                                   Optimization::VectorXf &gradient,
                                   Optimization::TripletListF &hessian)
            {
                twistEnergy(theta, energy);
                twistGradient(theta, gradient);
                return twistHessian(theta, hessian);
            };

        for (size_t step = 0; step < max_newton_iterations; ++step)
        {
            /** TODO: If we have different boundary conditions, we need to only pass the varying values here */
            Eigen::VectorXf theta = m_edge_theta;
            auto result = opt.step(theta);
            switch (result)
            {
            case Optimization::OptimizationStatus::SUCCESS:
                m_edge_theta = theta;
                break;
            case Optimization::OptimizationStatus::FAILURE:
                std::cout << "Material frame calculation failed" << std::endl;
                return;
            case Optimization::OptimizationStatus::CONVERGED:
                return;
            }
        }

        std::cout << "Reached max iterations in material frame optimization" << std::endl;
    }
}
