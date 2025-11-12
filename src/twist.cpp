#include <DiscreteElasticRod.hpp>
#include <Optimization.h>

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
        Optimization opt;
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
            auto theta = m_edge_theta;
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
    }
}
