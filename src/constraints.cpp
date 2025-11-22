#include <DiscreteElasticRod.hpp>
#include <Optimization.h>
#include <common.h>

#include <stdexcept>

bool DiscreteElasticRod::getConstraints(const Optimization::VectorXf &x_lambda, double &energy) const {
    print_debug("DiscreteElasticRod::getConstraints");

    //Eigen::VectorXf m_mass = Eigen::VectorXf::Ones(m_vertex_positions.size()) / 3.f;

    // Eigen::MatrixXf M = Eigen::MatrixXf::Zero(3 * m_n + 12, 3 * m_n + 12);
    // M.block<3, 3>(0, 0) = 4.f * Eigen::MatrixXf::Identity(3, 3);
    // M.block<3, 3>(3, 3) = m_mass.sum() * Eigen::MatrixXf::Identity(3, 3);

    // M.diagonal().tail(3 * m_n + 6) = m_mass;

    // Eigen::MatrixXf M = Eigen::MatrixXf::Zero(3 * (m_n + 2), 3 * (m_n + 2));
    // M.diagonal() = m_mass;

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

    // Eigen::VectorXf y = m_vertex_velocities.reshaped(3 * (m_n + 2), 1);

    Eigen::VectorXf C = Eigen::VectorXf::Zero(m_n + 1);
    for (int i = 0; i <= m_n; i++) {
        float new_edge_length = (x_lambda.segment<3>(3 * (i+1)) - x_lambda.segment<3>(3 * i)).norm();
        // $E_s = 1/2 (|x_{i+1} - x_i| - |\bar{e_i}|)^2$
        C(i) = 0.5f * std::pow(new_edge_length - m_edge_length(i), 2);
    }
    // C(m_n + 1) = q.squaredNorm() - 1.f;
    // //TODO: do we need to save rest positions for all vertices?
    // C.segment<3>(m_n + 2) = q._transformVector(m_vertex_positions.segment<3>(0)) + r - m_vertex_positions.segment<3>(0);
    // C.segment<3>(m_n + 5) = q._transformVector(m_vertex_positions.segment<3>(1)) + r - m_vertex_positions.segment<3>(1);

    //TODO: how do we get the values for lambda? is that passed in with the optimization vector?
    // Eigen::VectorXf lambda = x_lambda.tail(m_n + 1);

    //this is assuming that the output value is a scalar (using y^TMy instead of the yMy^T written in the paper), because otherwise the dimensions of this calculation don't work
    // energy = 0.5f * (x_lambda.head(3 * (m_n + 2)) - m_vertex_positions.reshaped(3 * (m_n + 2), 1)).squaredNorm() - C.dot(lambda);
    energy = C.sum();
    print_debug("the constraint energy is" + std::to_string(energy));
    return true;
}

bool DiscreteElasticRod::getConstraintGradient(const Optimization::VectorXf &x_lambda, Optimization::VectorXf &gradient) const {
    print_debug("DiscreteElasticRod::getConstraintGradient");

    auto maple_gradient = [](const Eigen::Vector3f &x_min_1, const Eigen::Vector3f &x_i, const Eigen::Vector3f &x_plus_1, const Eigen::Vector2f &e_bar)
    {
        Eigen::Vector3f dC_Stretch;
        dC_Stretch.setZero();

        float t1 = x_min_1[0] - x_i[0];
        float t2 = x_min_1[1] - x_i[1];
        float t3 = x_min_1[2] - x_i[2];
        float t4 = pow(t1, 0.2e1) + pow(t2, 0.2e1) + pow(t3, 0.2e1);
        float t5 = pow(t4, -0.1e1 / 0.2e1);
        float t6 = -x_plus_1[0] + x_i[0];
        float t7 = -x_plus_1[1] + x_i[1];
        float t8 = -x_plus_1[2] + x_i[2];
        float t9 = pow(t6, 0.2e1) + pow(t7, 0.2e1) + pow(t8, 0.2e1);
        float t10 = pow(t9, -0.1e1 / 0.2e1);
        t9 = (t9 * t10 - e_bar[1]) * t10;
        t4 = (t4 * t5 - e_bar[0]) * t5;
        dC_Stretch[0] = -t4 * t1 + t9 * t6;
        dC_Stretch[1] = -t4 * t2 + t9 * t7;
        dC_Stretch[2] = -t4 * t3 + t9 * t8;

        return dC_Stretch;
    };

    gradient = Eigen::VectorXf::Zero(3 * (m_n + 2));

    // For inextensibility, do not allow vertices 0 and n+1 to move, since they are fixed.
    // TODO: Validate if this sum range is correct
    for (int i = 1; i <= m_n; i++) {
        gradient.segment<3>(3 * i) = maple_gradient(
            x_lambda.segment<3>(3 * (i-1)),
            x_lambda.segment<3>(3 * i),
            x_lambda.segment<3>(3 * (i+1)),
            Eigen::Vector2f(m_edge_length(i-1), m_edge_length(i)));
    }

    print_debug("constraint gradient magnitude is " + std::to_string(gradient.norm()));
    return true;
}


bool DiscreteElasticRod::getConstraintHessian(const Optimization::VectorXf &x_lambda, Optimization::TripletListF &hessian) const {
    print_debug("DiscreteElasticRod::getConstraintHessian");

    auto maple_hessian = [](const Eigen::Vector3f &x_min_1, const Eigen::Vector3f &x_i, const Eigen::Vector3f &x_plus_1, const Eigen::Vector2f &e_bar)
    {
        Eigen::Matrix<float, 9, 1> d2C_Stretch;
        d2C_Stretch.setZero();

        float t1 = x_min_1[0] - x_i[0];
        float t2 = x_min_1[1] - x_i[1];
        float t3 = x_min_1[2] - x_i[2];
        float t4 = pow(t1, 0.2e1);
        float t5 = pow(t2, 0.2e1);
        float t6 = pow(t3, 0.2e1);
        float t7 = t6 + t5 + t4;
        float t8 = pow(t7, -0.3e1 / 0.2e1);
        float t9 = -pow(t7, 0.2e1) * t8 + e_bar[0];
        float t10 = t7 * t8;
        float t11 = x_i[0] - x_plus_1[0];
        float t12 = x_i[1] - x_plus_1[1];
        float t13 = x_i[2] - x_plus_1[2];
        float t14 = pow(t11, 0.2e1);
        float t15 = pow(t12, 0.2e1);
        float t16 = pow(t13, 0.2e1);
        float t17 = t16 + t15 + t14;
        float t18 = pow(t17, -0.3e1 / 0.2e1);
        float t19 = -pow(t17, 0.2e1) * t18 + e_bar[1];
        float t20 = t17 * t18;
        t7 = 0.1e1 / t7;
        t17 = 0.1e1 / t17;
        float t21 = t18 * t19 + t17;
        float t22 = t8 * t9 + t7;
        t2 = t22 * t2;
        float t23 = t12 * t11 * t21 + t2 * t1;
        t1 = t22 * t3 * t1 + t13 * t11 * t21;
        t2 = t13 * t12 * t21 + t2 * t3;
        d2C_Stretch[0] = t14 * t17 + t19 * (t14 * t18 - t20) + t4 * t7 + t9 * (t4 * t8 - t10);
        d2C_Stretch[1] = t23;
        d2C_Stretch[2] = t1;
        d2C_Stretch[3] = t23;
        d2C_Stretch[4] = t15 * t17 + t19 * (t15 * t18 - t20) + t5 * t7 + t9 * (t5 * t8 - t10);
        d2C_Stretch[5] = t2;
        d2C_Stretch[6] = t1;
        d2C_Stretch[7] = t2;
        d2C_Stretch[8] = t16 * t17 + t19 * (t16 * t18 - t20) + t6 * t7 + t9 * (t6 * t8 - t10);

        return d2C_Stretch.reshaped(3, 3);
    };

    for (int i = 1; i <= m_n; i++)
    {
        Eigen::Matrix3f xi_hessian = maple_hessian(
            x_lambda.segment<3>(3 * (i-1)),
            x_lambda.segment<3>(3 * i),
            x_lambda.segment<3>(3 * (i+1)),
            Eigen::Vector2f(m_edge_length(i-1), m_edge_length(i)));

        for (int offset = 0; offset < 3; offset++)
        {
            int j = i + offset - 1;
            for (int axis = 0; axis < 3; axis++)
            {
                hessian.emplace_back(i * 3 + axis, j * 3 + axis, xi_hessian(offset, axis));
            }
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
    print_debug(m_vertex_positions);

    print_debug("Before entering loop");
    for (size_t step = 0; step < max_newton_iterations; ++step) {
        print_debug("At start of loop" + std::to_string(step));
        Eigen::VectorXf x_lambda = Eigen::VectorXf::Zero(3 * (m_n + 2) + m_n + 1);
        print_debug("X_lambda created");
        x_lambda.head(3 * (m_n + 2)) = m_vertex_positions.reshaped(3 * (m_n + 2), 1);
        print_debug("X_lambda head set");
        x_lambda.tail(m_n + 1) = Eigen::VectorXf::Constant(m_n + 1, 1.f);
        print_debug("X_lambda tail set");

        print_debug("After initialization in loop " + std::to_string(step));
        auto result = opt.step(x_lambda);
        switch (result)
        {
            case Optimization::OptimizationStatus::SUCCESS:
                print_debug("Start set vertex positions");
                m_vertex_positions = x_lambda.head(3 * (m_n + 2)).reshaped(3, m_n + 2);
                print_debug("End set vertex positions to ");
                print_debug(m_vertex_positions);
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
