#include <DiscreteElasticRod.hpp>
#include <Optimization.h>
#include <common.h>

#include <stdexcept>

#if false
    #define FIRST_METHOD
#else
    #define SECOND_METHOD
#endif

#if false
    #define ENDPOINTS
#endif

bool DiscreteElasticRod::getConstraints(const Opt::VectorX &x_lambda, Opt::Float &energy) const {
    print_debug("");

    VectorX C = VectorX::Zero(m_n + 1);
    for (int i = 0; i <= m_n; i++) {
        //For first constraint formulation
        // $E_s = 1/2 (|x_{i+1} - x_i| - |\bar{e_i}|)^2$
#ifdef FIRST_METHOD
        Float new_edge_length = (x_lambda.segment<3>(3 * (i+1)) - x_lambda.segment<3>(3 * i)).norm();
        C(i) = 0.5f * std::pow(new_edge_length - m_edge_length(i), 2);
#endif

#ifdef SECOND_METHOD
        //For second constraint formulation
        // $E_s = 1/2 (|x_{i+1} - x_i|^2 - |\bar{e_i}|^2)$
        Float squared_edge_length = (x_lambda.segment<3>(3 * (i + 1)) - x_lambda.segment<3>(3 * i)).squaredNorm();
        C(i) = 0.5 * (squared_edge_length - std::pow(m_edge_length(i), 2));
        if (C(i) > 100) {
            print_debug("vertex " + std::to_string(i+1) + ": ");
            print_debug(x_lambda.segment<3>(3 * (i + 1)));
            print_debug("vertex " + std::to_string(i) + ": ");
            print_debug(x_lambda.segment<3>(3 * i));
            print_debug("because squared_edge_length is " + std::to_string(squared_edge_length) + " and m_edge_length is " + std::to_string(m_edge_length(i)));
            print_debug("energy at vertex " + std::to_string(i) + " is too high" );
        }
#endif

    }

    energy = C.sum();

    if (energy > 100) {
        print_debug("energy too high wtf");
    }
    print_debug("the constraint energy is " + std::to_string(energy));
    return true;
}

bool DiscreteElasticRod::getConstraintGradient(const Opt::VectorX &x_lambda, Opt::VectorX &gradient) const {
    print_debug("");

    auto maple_gradient = [](const Vector3 &x_min_1, const Vector3 &x_i, const Vector3 &x_plus_1, const Vector2 &e_bar)
    {
        Vector3 dC_Stretch;
        dC_Stretch.setZero();

#ifdef FIRST_METHOD
        //For first constraint formulation
        Float t1 = x_min_1[0] - x_i[0];
        Float t2 = x_min_1[1] - x_i[1];
        Float t3 = x_min_1[2] - x_i[2];
        Float t4 = pow(t1, 0.2e1) + pow(t2, 0.2e1) + pow(t3, 0.2e1);
        Float t5 = pow(t4, -0.1e1 / 0.2e1);
        Float t6 = -x_plus_1[0] + x_i[0];
        Float t7 = -x_plus_1[1] + x_i[1];
        Float t8 = -x_plus_1[2] + x_i[2];
        Float t9 = pow(t6, 0.2e1) + pow(t7, 0.2e1) + pow(t8, 0.2e1);
        Float t10 = pow(t9, -0.1e1 / 0.2e1);
        t9 = (t9 * t10 - e_bar[1]) * t10;
        t4 = (t4 * t5 - e_bar[0]) * t5;
        dC_Stretch[0] = -t4 * t1 + t9 * t6;
        dC_Stretch[1] = -t4 * t2 + t9 * t7;
        dC_Stretch[2] = -t4 * t3 + t9 * t8;
#endif


#ifdef SECOND_METHOD
        //For second constraint formulation
        Float t1 = 2;
        dC_Stretch[0] = t1 * x_i[0] - x_min_1[0] - x_plus_1[0];
        dC_Stretch[1] = t1 * x_i[1] - x_min_1[1] - x_plus_1[1];
        dC_Stretch[2] = t1 * x_i[2] - x_min_1[2] - x_plus_1[2];
#endif

        return dC_Stretch;
    };

#ifdef ENDPOINTS
    //for first vertex
    auto maple_first_vertex_gradient = [](const Vector3 &x_plus_1, const Vector3 &x_i, const Vector2 &e_bar) {
        Vector3 dC_Stretch_plus;
        dC_Stretch_plus.setZero();

        // For first constraint formulation
#ifdef FIRST_METHOD
        Float t1 = x_plus_1[0] - x_i[0];
        Float t2 = x_plus_1[1] - x_i[1];
        Float t3 = x_plus_1[2] - x_i[2];
        Float t4 = pow(t1, 0.2e1) + pow(t2, 0.2e1) + pow(t3, 0.2e1);
        Float t5 = pow(t4, -0.1e1 / 0.2e1);
        t4 = (-t4 * t5 + e_bar[1]) * t5;
        dC_Stretch_plus[0] = t4 * t1;
        dC_Stretch_plus[1] = t4 * t2;
        dC_Stretch_plus[2] = t4 * t3;
#endif

#ifdef SECOND_METHOD
        //For second constraint formulation
        dC_Stretch_plus[0] = -x_plus_1[0] + x_i[0];
        dC_Stretch_plus[1] = -x_plus_1[1] + x_i[1];
        dC_Stretch_plus[2] = -x_plus_1[2] + x_i[2];
#endif

        return dC_Stretch_plus;
    };

    //For last vertex
    auto maple_last_vertex_gradient = [](const Vector3 &x_min_1, const Vector3 &x_i, const Vector2 &e_bar) {
        Vector3 dC_Stretch_min;
        dC_Stretch_min.setZero();

#ifdef FIRST_METHOD
        // For first constraint formulation
        Float t1 = x_min_1[0] - x_i[0];
        Float t2 = x_min_1[1] - x_i[1];
        Float t3 = x_min_1[2] - x_i[2];
        Float t4 = pow(t1, 0.2e1) + pow(t2, 0.2e1) + pow(t3, 0.2e1);
        Float t5 = pow(t4, -0.1e1 / 0.2e1);
        t4 = (t4 * t5 - e_bar[0]) * t5;
        dC_Stretch_min[0] = -t4 * t1;
        dC_Stretch_min[1] = -t4 * t2;
        dC_Stretch_min[2] = -t4 * t3;
#endif

#ifdef SECOND_METHOD
        //For second constraint formulation
        dC_Stretch_min[0] = x_i[0] - x_min_1[0];
        dC_Stretch_min[1] = x_i[1] - x_min_1[1];
        dC_Stretch_min[2] = x_i[2] - x_min_1[2];
#endif

        return dC_Stretch_min;
    };
#endif

    gradient = VectorX::Zero(3 * (m_n + 2));

    // For inextensibility, do not allow vertices 0 and n+1 to move, since they are fixed.
    // TODO: Validate if this sum range is correct
    for (int i = 1; i <= m_n; i++) {
        gradient.segment<3>(3 * i) = maple_gradient(
            x_lambda.segment<3>(3 * (i-1)),
            x_lambda.segment<3>(3 * i),
            x_lambda.segment<3>(3 * (i+1)),
            Vector2(m_edge_length(i-1), m_edge_length(i)));
    }

#ifdef ENDPOINTS
    //For first and last vertices
     gradient.segment<3>(0) = maple_first_vertex_gradient(
         x_lambda.segment<3>(3),
             x_lambda.segment<3>(0),
             Vector2(0.f, m_edge_length(0)));
     gradient.segment<3>(3 * (m_n + 1)) = maple_last_vertex_gradient(
         x_lambda.segment<3>(3 * m_n),
             x_lambda.segment<3>(3 * (m_n + 1)),
                 Vector2(m_edge_length(m_n), 0.f));
#endif


    print_debug("constraint gradient magnitude is " + std::to_string(gradient.norm()));
    return true;
}


bool DiscreteElasticRod::getConstraintHessian(const Opt::VectorX &x_lambda, Opt::TripletList &hessian) const {
    print_debug("DiscreteElasticRod::getConstraintHessian");

    auto maple_hessian = [](const Vector3 &x_min_1, const Vector3 &x_i, const Vector3 &x_plus_1, const Vector2 &e_bar)
    {
        Eigen::Matrix<Float, 9, 1> d2C_Stretch;
        d2C_Stretch.setZero();

        //For first constraint formulation
#ifdef FIRST_METHOD
        Float t1 = x_min_1[0] - x_i[0];
        Float t2 = x_min_1[1] - x_i[1];
        Float t3 = x_min_1[2] - x_i[2];
        Float t4 = pow(t1, 0.2e1);
        Float t5 = pow(t2, 0.2e1);
        Float t6 = pow(t3, 0.2e1);
        Float t7 = t6 + t5 + t4;
        Float t8 = pow(t7, -0.3e1 / 0.2e1);
        Float t9 = -pow(t7, 0.2e1) * t8 + e_bar[0];
        Float t10 = t7 * t8;
        Float t11 = x_i[0] - x_plus_1[0];
        Float t12 = x_i[1] - x_plus_1[1];
        Float t13 = x_i[2] - x_plus_1[2];
        Float t14 = pow(t11, 0.2e1);
        Float t15 = pow(t12, 0.2e1);
        Float t16 = pow(t13, 0.2e1);
        Float t17 = t16 + t15 + t14;
        Float t18 = pow(t17, -0.3e1 / 0.2e1);
        Float t19 = -pow(t17, 0.2e1) * t18 + e_bar[1];
        Float t20 = t17 * t18;
        t7 = 0.1e1 / t7;
        t17 = 0.1e1 / t17;
        Float t21 = t18 * t19 + t17;
        Float t22 = t8 * t9 + t7;
        t2 = t22 * t2;
        Float t23 = t12 * t11 * t21 + t2 * t1;
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
#endif

        //For second constraint formulation
#ifdef SECOND_METHOD
         d2C_Stretch[0] = 2;
         d2C_Stretch[1] = 0;
         d2C_Stretch[2] = 0;
         d2C_Stretch[3] = 0;
         d2C_Stretch[4] = 2;
         d2C_Stretch[5] = 0;
         d2C_Stretch[6] = 0;
         d2C_Stretch[7] = 0;
         d2C_Stretch[8] = 2;
#endif

        return d2C_Stretch.reshaped(3, 3);
    };

#ifdef ENDPOINTS
    //For first vertex
    auto maple_hessian_last_vertex = [](const Vector3 &x_plus_1, const Vector3 &x_i, const Vector2 &e_bar) {
        Eigen::Matrix<Float, 9, 1> d2C_Stretch_plus;
        d2C_Stretch_plus.setZero();

        //For first constraint formulation
        Float t1 = x_plus_1[0] - x_i[0];
        Float t2 = x_plus_1[1] - x_i[1];
        Float t3 = x_plus_1[2] - x_i[2];
        Float t4 = pow(t2, 0.2e1);
        Float t5 = pow(t1, 0.2e1);
        Float t6 = pow(t3, 0.2e1);
        Float t7 = t4 + t5 + t6;
        Float t8 = pow(t7, -0.3e1 / 0.2e1);
        Float t9 = -pow(t7, 0.2e1) * t8 + e_bar[1];
        Float t10 = t7 * t8;
        t7 = 0.1e1 / t7;
        Float t11 = t8 * t9 + t7;
        Float t12 = t2 * t1 * t11;
        t1 = t3 * t1 * t11;
        t2 = t3 * t2 * t11;
        d2C_Stretch_plus[0] = t5 * t7 + t9 * (t5 * t8 - t10);
        d2C_Stretch_plus[1] = t12;
        d2C_Stretch_plus[2] = t1;
        d2C_Stretch_plus[3] = t12;
        d2C_Stretch_plus[4] = t4 * t7 + t9 * (t4 * t8 - t10);
        d2C_Stretch_plus[5] = t2;
        d2C_Stretch_plus[6] = t1;
        d2C_Stretch_plus[7] = t2;
        d2C_Stretch_plus[8] = t6 * t7 + t9 * (t6 * t8 - t10);

        //In second constraint formulation, just Identity

        return d2C_Stretch_plus.reshaped(3, 3);
    };

    auto maple_hessian_first_vertex = [](const Vector3 &x_min_1, const Vector3 &x_i, const Vector2 &e_bar) {
        Eigen::Matrix<Float, 9, 1> d2C_Stretch_min;
        d2C_Stretch_min.setZero();

        //For first constraint formulation
        Float t1 = x_min_1[0] - x_i[0];
        Float t2 = x_min_1[1] - x_i[1];
        Float t3 = x_min_1[2] - x_i[2];
        Float t4 = pow(t3, 0.2e1);
        Float t5 = pow(t1, 0.2e1);
        Float t6 = pow(t2, 0.2e1);
        Float t7 = t4 + t5 + t6;
        Float t8 = pow(t7, -0.3e1 / 0.2e1);
        Float t9 = -pow(t7, 0.2e1) * t8 + e_bar[0];
        Float t10 = t7 * t8;
        t7 = 0.1e1 / t7;
        Float t11 = t8 * t9 + t7;
        Float t12 = t2 * t1 * t11;
        t1 = t3 * t1 * t11;
        t2 = t3 * t2 * t11;
        d2C_Stretch_min[0] = t5 * t7 + t9 * (t5 * t8 - t10);
        d2C_Stretch_min[1] = t12;
        d2C_Stretch_min[2] = t1;
        d2C_Stretch_min[3] = t12;
        d2C_Stretch_min[4] = t6 * t7 + t9 * (t6 * t8 - t10);
        d2C_Stretch_min[5] = t2;
        d2C_Stretch_min[6] = t1;
        d2C_Stretch_min[7] = t2;
        d2C_Stretch_min[8] = t4 * t7 + t9 * (t4 * t8 - t10);

        //For second constraint formulation, just identity

        return d2C_Stretch_min.reshaped(3, 3);
    };
#endif


    for (int i = 1; i <= m_n; i++)
    {
        Vector3 x_min_1 = x_lambda.segment<3>(3 * (i-1));
        Vector3 x_i = x_lambda.segment<3>(3 * i);
        Vector3 x_plus_1 = x_lambda.segment<3>(3 * (i+1));
        Vector2 e_bar = Vector2(m_edge_length(i-1), m_edge_length(i));

        Matrix3 xi_hessian = maple_hessian(
            x_min_1,
            x_i,
            x_plus_1,
            e_bar);

        for (int offset = 0; offset < 3; offset++)
        {
            int j = i + offset - 1;
            for (int axis = 0; axis < 3; axis++)
            {
                hessian.emplace_back(i * 3 + axis, j * 3 + axis, xi_hessian(offset, axis));
            }
        }
    }

#ifdef ENDPOINTS
    Matrix3 maple_hessian_first =
        //For first constraint formulation
#ifdef  FIRST_METHOD
        maple_hessian_first_vertex(
         x_lambda.segment<3>(3),
         x_lambda.segment<3>(0),
         Vector2(m_edge_length(0), 0.f));
#endif FIRST_METHOD

        //For second constraint formulation
#ifdef SECOND_METHOD
        Matrix3::Identity();
#endif SECOND_METHOD

    Matrix3 maple_hessian_last =
        //For first constraint formulation
#ifdef  FIRST_METHOD
        maple_hessian_last_vertex(
     x_lambda.segment<3>(3 * m_n),
        x_lambda.segment<3>(3 * (m_n + 1)),
        Vector2(0.f, m_edge_length(m_n)));
#endif

        //For second constraint formulation
#ifdef SECOND_METHOD
        Matrix3::Identity();
#endif

    // For endpoint constraints
    for (int offset = 0; offset < 3; offset++) {
        for (int axis = 0; axis < 3; axis++) {
            hessian.emplace_back(axis, axis, maple_hessian_first(offset, axis));
            hessian.emplace_back(axis, 3 + axis, maple_hessian_first(offset, axis));

            hessian.emplace_back(3 * m_n + axis,  3 * m_n + axis, maple_hessian_last(offset, axis));
            hessian.emplace_back(3 * m_n + axis,  3 * m_n + 3 + axis, maple_hessian_last(offset, axis));
        }
    }
#endif


    print_debug("[getHessian] exit");

    return true;
}

void DiscreteElasticRod::applyConstraints(size_t max_newton_iterations) {
    print_debug("[applyConstraints] enter");

    Opt opt;
    opt.optimizer = Opt::Optimizer::NEWTON;
    opt.tolerance_exponent = -10;
    opt.objective_function = [&](const Opt::VectorX &x_lambda, Opt::Float &energy)
        {
            return getConstraints(x_lambda, energy);
        };
    opt.gradient_function = [&](const Opt::VectorX &x_lambda,
                                Opt::Float &energy,
                                Opt::VectorX &gradient)
        {
            getConstraints(x_lambda, energy);
            return getConstraintGradient(x_lambda, gradient);
        };
    opt.hessian_function = [&](const Opt::VectorX &x_lambda,
                                Opt::Float &energy,
                                Opt::VectorX &gradient,
                                Opt::TripletList &hessian)
        {
            getConstraints(x_lambda, energy);
            getConstraintGradient(x_lambda, gradient);
            return getConstraintHessian(x_lambda, hessian);
        };

    print_debug("Initial positions:");
    print_debug(m_vertex_positions);
    print_debug("Initial edge length:");
    print_debug(VectorX(getEdges().colwise().norm()).transpose());
    print_debug("Initial edge length difference:");
    print_debug((VectorX(getEdges().colwise().norm()) - m_edge_length).transpose());

    bool exit = false;
    for (size_t step = 0; step < max_newton_iterations; ++step) {
        VectorX x_lambda = VectorX::Zero(3 * (m_n + 2));
        x_lambda.head(3 * (m_n + 2)) = m_vertex_positions.reshaped(3 * (m_n + 2), 1);

        auto result = opt.step(x_lambda);
        switch (result)
        {
            case Opt::OptimizationStatus::SUCCESS:
                m_vertex_positions = x_lambda.head(3 * (m_n + 2)).reshaped(3, m_n + 2);
                break;
            case Opt::OptimizationStatus::FAILURE:
                std::cout << "constraint calculation fail" << std::endl;
                exit = true;
                break;
            case Opt::OptimizationStatus::CONVERGED:
                exit = true;
                break;
        }

        if (exit)
        {
            break;
        }
    }

    print_debug("Resulting positions:");
    print_debug(m_vertex_positions);
    print_debug("Result edge length:");
    print_debug(VectorX(getEdges().colwise().norm()).transpose());
    print_debug("Result edge length difference:");
    print_debug((VectorX(getEdges().colwise().norm()) - m_edge_length).transpose());
    print_debug("exit");
}
