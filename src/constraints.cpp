#include <DiscreteElasticRod.hpp>
#include <Optimization.h>
#include <common.h>

#include <stdexcept>

bool DiscreteElasticRod::getConstraints(const Optimization::VectorXf &x_lambda, double &energy) const {
    print_debug("DiscreteElasticRod::getConstraints");

    Eigen::VectorXf m_mass = Eigen::VectorXf::Ones(m_vertex_positions.size()) / 3.f;

    // Eigen::MatrixXf M = Eigen::MatrixXf::Zero(3 * m_n + 12, 3 * m_n + 12);
    // M.block<3, 3>(0, 0) = 4.f * Eigen::MatrixXf::Identity(3, 3);
    // M.block<3, 3>(3, 3) = m_mass.sum() * Eigen::MatrixXf::Identity(3, 3);

    // M.diagonal().tail(3 * m_n + 6) = m_mass;

    Eigen::MatrixXf M = Eigen::MatrixXf::Zero(3 * (m_n + 2), 3 * (m_n + 2));
    M.diagonal() = m_mass;

    // Eigen::Quaternionf q = Eigen::Quaternionf(1, 0, 0, 0);

    // //TODO: get q_deriv apparently there is also a q.derived() function? but that can't do what i think it does, can it?
    // Eigen::Quaternionf q_deriv = Eigen::Quaternionf(1, 0, 0, 0);

    // Eigen::Vector3f r = Eigen::Vector3f::Zero();
    // //TODO: and how are we supposed to get the derivative of r in time either?
    // Eigen::Vector3f r_deriv = Eigen::Vector3f(1, 0, 0);

    // Eigen::VectorXf y = Eigen::VectorXf::Zero(3 * m_n + 12);
    // y.segment<3>(0) = (q.inverse() * q_deriv).coeffsScalarFirst().segment<3>(1);
    // y.segment<3>(3) = r_deriv;
    // y.tail(3 * m_n + 6) = m_vertex_velocities;

    Eigen::VectorXf y = m_vertex_velocities.reshaped(3 * (m_n + 2), 1);

    Eigen::VectorXf C = Eigen::VectorXf::Zero(m_n + 1);
    for (int i = 0; i <= m_n; i++) {
        C(i) = (x_lambda.segment<3>(3 * (i+1)) - x_lambda.segment<3>(3 * i)).squaredNorm() - m_edge_length(i) * m_edge_length(i);
    }
    // C(m_n + 1) = q.squaredNorm() - 1.f;
    // //TODO: do we need to save rest positions for all vertices?
    // C.segment<3>(m_n + 2) = q._transformVector(m_vertex_positions.segment<3>(0)) + r - m_vertex_positions.segment<3>(0);
    // C.segment<3>(m_n + 5) = q._transformVector(m_vertex_positions.segment<3>(1)) + r - m_vertex_positions.segment<3>(1);

    //TODO: how do we get the values for lambda? is that passed in with the optimization vector?
    Eigen::VectorXf lambda = x_lambda.tail(m_n + 1);

    //this is assuming that the output value is a scalar (using y^TMy instead of the yMy^T written in the paper), because otherwise the dimensions of this calculation don't work
    energy = 0.5f * (x_lambda.head(3 * (m_n + 2)) - m_vertex_positions).squaredNorm() - C.dot(lambda); // TODO: Cannot take squared norm on matrix

    return true;
}

bool DiscreteElasticRod::getConstraintGradient(const Optimization::VectorXf &x_lambda, Optimization::VectorXf &gradient) const {
    print_debug("DiscreteElasticRod::getConstraintGradient");

    Eigen::VectorXf X_constraints = Eigen::VectorXf::Zero(3 * (m_n + 2));

    //TODO: Figure out how to take the derivative of y^TMy - Cl

    // contribution from loss function
    Eigen::VectorXf lambdas = x_lambda.tail(m_n + 1);
    X_constraints = x_lambda.head(3 * (m_n + 2)) - m_vertex_positions;

    // contribution from inextensibility constraints
    for (int i = 0; i <= m_n; i++) {
        X_constraints.segment<3>(3 * i) += 2 * lambdas(i) * (x_lambda.segment<3>(3 * i) - x_lambda.segment<3>(3 * (i + 1)));
    }
    for (int i = 1; i <= m_n + 1; i++) {
        X_constraints.segment<3>(3 * i) += 2 * lambdas(i-1) * (x_lambda.segment<3>(3 * i) - x_lambda.segment<3>(3 * (i - 1)));
    }

    Eigen::VectorXf C = Eigen::VectorXf::Zero(m_n + 1);
    for (int i = 0; i <= m_n; i++) {
        C(i) = (x_lambda.segment<3>(3 * (i+1)) - x_lambda.segment<3>(3 * i)).squaredNorm() - m_edge_length(i) * m_edge_length(i);
    }
    gradient = Eigen::VectorXf::Zero(3 * (m_n + 2) + m_n + 1);

    gradient.head(3 * (m_n + 2)) = X_constraints;
    gradient.tail(m_n + 1) = C;

    return true;
}


bool DiscreteElasticRod::getConstraintHessian(const Optimization::VectorXf &x_lambda, Optimization::TripletListF &hessian) const {
    print_debug("DiscreteElasticRod::getConstraintHessian");

    Eigen::VectorXf x = x_lambda.head(3 * (m_n + 2));
    Eigen::VectorXf lambda = x_lambda.tail(m_n + 1);


    for (int i = 0; i <= m_n+1; i++) {
        for (int j = 0; j < 3; j++) {
            if (i == 0) {
                hessian.emplace_back(0 + j, 0 + j, 1.f + 2.f * lambda(0));
                hessian.emplace_back(0 + j, 3 + j, -2.f * lambda(0));

            } else if (i == m_n + 1) {
                hessian.emplace_back(3 * (m_n + 1) + j, 3 * (m_n + 1) + j, 1.f + 2.f * lambda(m_n));
                hessian.emplace_back(3 * (m_n + 1) + j, 3 * m_n + j, -2.f * lambda(m_n));

            } else {
                hessian.emplace_back(3 * i + j, 3 * i + j, 1.f + 2.f * lambda(i) + 2.f * lambda(i-1));
                hessian.emplace_back(3 * i + j, 3 * (i+1) + j, - 2.f * lambda(i));
                hessian.emplace_back(3 * i + j, 3 * (i-1) + j, - 2.f * lambda(i-1));
            }
        }
    }

    for (int i = 0; i <= m_n; i++) {
        for (int j = 0; j < 3; j++) {
            int dl_index = 3 * (m_n + 2 + i) + j;
            int dx_index = 3 * i + j;

            float x_i_minus_x_i_plus_one = 2 * (x(3 * i + j) - x(3 * (i + 1) + j));
            float x_i_plus_one_minus_x_i = 2 * (x(3 * (i + 1) + j) - x(3 * i + j));

            hessian.emplace_back(dl_index, dx_index, x_i_minus_x_i_plus_one);
            hessian.emplace_back(dl_index, dx_index + 3, x_i_plus_one_minus_x_i);

            hessian.emplace_back(dx_index, dl_index, x_i_minus_x_i_plus_one);
            hessian.emplace_back(dx_index + 3, dl_index, x_i_plus_one_minus_x_i);
        }

    }

    return true;
}

void DiscreteElasticRod::applyConstraints(size_t max_newton_iterations) {
    Optimization opt;
    opt.optimizer = Optimization::Optimizer::NEWTON;
    opt.tolerance_exponent = -10;
    opt.objective_function = [&](const Eigen::VectorXf &x_lambda, double &energy)
        {
            return getConstraints(x_lambda, energy);
        };
    opt.gradient_function = [&](const Eigen::VectorXf &x_lambda,
                                double &energy,
                                Eigen::VectorXf &gradient)
        {
            getConstraints(x_lambda, energy);
            return getConstraintGradient(x_lambda, gradient);
        };
    opt.hessian_function = [&](const Optimization::VectorXf &x_lambda,
                                double &energy,
                                Optimization::VectorXf &gradient,
                                Optimization::TripletListF &hessian)
        {
            getConstraints(x_lambda, energy);
            getConstraintGradient(x_lambda, gradient);
            return getConstraintHessian(x_lambda, hessian);
        };

    std::cout << "before entering loop" << std::endl;

    print_debug("Before entering loop");
    for (size_t step = 0; step < max_newton_iterations; ++step) {
        print_debug("At start of loop" + std::to_string(step));
        Eigen::VectorXf x_lambda = Eigen::VectorXf::Zero(3 * (m_n + 2) + m_n + 1);
        print_debug("X_lambda created");
        x_lambda.head(3 * (m_n + 2)) = m_vertex_positions.reshaped(3 * (m_n + 2), 1);
        print_debug("X_lambda head set");
        x_lambda.tail(m_n + 1) = Eigen::VectorXf::Constant(m_n + 1, 0.5f);
        print_debug("X_lambda tail set");

        print_debug("After initialization in loop " + std::to_string(step));
        auto result = opt.step(x_lambda); // TODO: This crashes
        switch (result)
        {
            case Optimization::OptimizationStatus::SUCCESS:
                print_debug("Start set vertex positions");
                m_vertex_positions = x_lambda.head(3 * (m_n + 2)).reshaped(3, m_n + 2);
                print_debug("End set vertex positions");
                break;
            case Optimization::OptimizationStatus::FAILURE:
                std::cout << "constraint calculation fail" << std::endl;
                return;
            case Optimization::OptimizationStatus::CONVERGED:
                return;
        }
    }

    std::cout << "Reached max iterations in constraint calculation" << std::endl;

}
