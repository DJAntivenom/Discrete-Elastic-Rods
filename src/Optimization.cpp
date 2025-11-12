#include <Optimization.h>

#include <iostream>

Optimization::OptimizationStatus Optimization::step(VectorXf &y)
{

    VectorXf dy;
    double initial_objective_value;
    OptimizationStatus status;
    switch (optimizer)
    {
    case Optimization::Optimizer::GRADIENT_DESCENT:
        status = getDirectionGradientDescent(y, dy, initial_objective_value);
        break;
    case Optimization::Optimizer::NEWTON:
        status = getDirectionNewton(y, dy, initial_objective_value);
        break;
    default:
        assert(0);
    }

    if (status == OptimizationStatus::SUCCESS)
    {
        status = lineSearch(y, dy, initial_objective_value) ? OptimizationStatus::SUCCESS : OptimizationStatus::FAILURE;
    }

    return status;
}

Optimization::OptimizationStatus Optimization::getDirectionGradientDescent(const VectorXf &y, VectorXf &dy,
                                                                           double &initial_objective_value)
{
    VectorXf gradient;
    if (!gradient_function(y, initial_objective_value, gradient))
    {
        assert(0);
        return OptimizationStatus::FAILURE;
    }

    if (gradient.squaredNorm() / y.rows() < pow(10.0, tolerance_exponent))
        return OptimizationStatus::CONVERGED;

    dy = -gradient;
    return OptimizationStatus::SUCCESS;
}

Optimization::OptimizationStatus Optimization::getDirectionNewton(const VectorXf &y, VectorXf &dy,
                                                                  double &initial_objective_value)
{
    VectorXf gradient;
    TripletListF hessian;

    if (!hessian_function(y, initial_objective_value, gradient, hessian))
    {
        assert(0);
        return OptimizationStatus::FAILURE;
    }

    if (gradient.squaredNorm() / y.rows() < pow(10.0, tolerance_exponent))
        return OptimizationStatus::CONVERGED;

    /// Linear solve to obtain search direction.
    bool linear_solve_success;
    linear_solve_success = linearSolve(hessian, -gradient, dy);

    return linear_solve_success ? OptimizationStatus::SUCCESS : OptimizationStatus::FAILURE;
}

bool Optimization::lineSearch(VectorXf &y, const VectorXf &dy, double initial_objective_value)
{
    double alpha = 1.0;

    int max_steps = 50;
    for (int i = 0; i < max_steps; i++)
    {
        VectorXf y_line_search = y + alpha * dy;

        double new_objective_value;
        if (objective_function(y_line_search, new_objective_value))
        {
            if (new_objective_value < initial_objective_value)
            {
                y = y_line_search;
                return true;
            }
        }
        else
        {
            /// Invalid y. Go to next loop iteration.
        }
        alpha *= 0.5;
    }

    return false;
}

bool Optimization::lineSearch(VectorXf &y, const VectorXf &dy)
{
    double initial_objective_value;
    objective_function(y, initial_objective_value);
    return lineSearch(y, dy, initial_objective_value);
}

bool Optimization::linearSolve(const TripletListF &hessian, const VectorXf &b, VectorXf &x)
{
    Solver solver;

    SparseMatrixF A(b.rows(), b.rows());
    A.setFromTriplets(hessian.begin(), hessian.end());

    /// Start with very small diagonal regularizer to ensure diagonal is included in sparsity pattern.
    SparseMatrixF H(A.rows(), A.cols());
    H.setIdentity();
    H.diagonal().array() = 1e-10;
    solver.analyzePattern(A + H);

    int indefinite_count_reg_cnt = 0, invalid_search_dir_cnt = 0, invalid_residual_cnt = 0;
    double alpha = 1e-6;
    for (int i = 0; i < 50; i++)
    {
        solver.factorize(A + H);

        /// Should have NumericalIssue if K is not positive definite.
        if (solver.info() == Eigen::NumericalIssue)
        {
            H.diagonal().array() = alpha;
            alpha *= 10;
            indefinite_count_reg_cnt++;
            continue;
        }

        x = solver.solve(b);

        /// Check that step has positive dot product with residual.
        double dot_x_b = x.normalized().dot(b.normalized());
        bool search_dir_correct_sign = dot_x_b > 0;
        if (!search_dir_correct_sign)
        {
            invalid_search_dir_cnt++;
        }

        /// Check for reasonable step size.
        bool solve_success = x.norm() < 1e3;
        if (!solve_success)
        {
            invalid_residual_cnt++;
        }

        /// In case of problem, increase regularizer by an order of magnitude and try again.
        if (search_dir_correct_sign && solve_success)
        {
            return true;
        }
        else
        {
            H.diagonal().array() = alpha;
            alpha *= 10;
        }
    }

    std::cout << "Linear solve failure reasons: " << indefinite_count_reg_cnt << " indefinite, "
        << invalid_search_dir_cnt << " invalid search dir, " << invalid_residual_cnt << " invalid residual."
        << std::endl;

    return false;
}

bool Optimization::checkGradient(const VectorXf &y, VectorXf &error,
                                 const std::function<bool(const VectorXf &, double &)> &func,
                                 const std::function<bool(const VectorXf &, VectorXf &)> &grad_func, double epsilon,
                                 int print_level)
{
    int num_variables = (int)y.rows();

    VectorXf grad;
    if (!(grad_func(y, grad)))
    {
        /// Gradient evaluation failed.
        error = VectorXf::Constant(num_variables, NAN);
        if (print_level > 0)
            std::cout << "Optimization::checkGradient: Gradient evaluation failed at given state!" << std::endl;
        return false;
    }

    VectorXf grad_fd(num_variables);
    for (int i = 0; i < num_variables; ++i)
    {
        VectorXf y_plus = y, y_minus = y, y_plus2 = y, y_minus2 = y;
        y_plus(i) += epsilon;
        y_minus(i) -= epsilon;
        y_plus2(i) += 2.0 * epsilon;
        y_minus2(i) -= 2.0 * epsilon;

        double f_plus, f_minus, f_plus2, f_minus2;
        if (func(y_plus, f_plus) && func(y_minus, f_minus) && func(y_plus2, f_plus2) && func(y_minus2, f_minus2))
        {
            double grad_fd_i = (f_plus - f_minus) / (2.0 * epsilon);
            double error_ratio =
                ((f_plus2 - f_minus2) / (4.0 * epsilon) - grad(i)) / ((f_plus - f_minus) / (2.0 * epsilon) - grad(i));

            if ((print_level == 1 && fabs(grad_fd_i - grad(i)) > epsilon) || print_level > 1)
                std::cout << "Optimization::checkGradient: grad_fd[" << i << "] = " << grad_fd_i << ", grad[" << i
                << "] = " << grad(i) << ", error ratio = " << error_ratio << std::endl;

            grad_fd(i) = grad_fd_i;
        }
        else
        {
            if (print_level > 0)
                std::cout << "Optimization::checkGradient: Function evaluation failed!" << std::endl;
            grad_fd(i) = NAN;
        }
    }

    error = grad - grad_fd;
    return error.allFinite() && error.cwiseAbs().maxCoeff() < epsilon;
}

bool Optimization::checkHessian(const VectorXf &y, MatrixXf &error,
                                const std::function<bool(const VectorXf &, VectorXf &)> &grad_func,
                                const std::function<bool(const VectorXf &, MatrixXf &)> &hess_func, double epsilon,
                                int print_level)
{
    int num_variables = (int)y.rows();

    MatrixXf hess;
    if (!(hess_func(y, hess)))
    {
        /// Hessian function failed.
        if (print_level > 0)
            std::cout << "Optimization::checkHessian: Hessian evaluation failed at given state!" << std::endl;
        error = MatrixXf::Constant(num_variables, num_variables, NAN);
        return false;
    }

    MatrixXf hess_fd(num_variables, num_variables);
    for (int i = 0; i < num_variables; ++i)
    {
        VectorXf y_plus = y, y_minus = y, y_plus2 = y, y_minus2 = y;
        y_plus(i) += epsilon;
        y_minus(i) -= epsilon;
        y_plus2(i) += 2.0 * epsilon;
        y_minus2(i) -= 2.0 * epsilon;

        VectorXf grad_plus, grad_minus, grad_plus2, grad_minus2;
        if (grad_func(y_plus, grad_plus) && grad_func(y_minus, grad_minus) && grad_func(y_plus2, grad_plus2) &&
            grad_func(y_minus2, grad_minus2))
        {
            for (int j = i; j < num_variables; ++j)
            { /// Upper triangular hessian
                double hess_fd_i_j = (grad_plus(j) - grad_minus(j)) / (2.0 * epsilon);
                double error_ratio = ((grad_plus2(j) - grad_minus2(j)) / (4.0 * epsilon) - hess(i, j)) /
                    ((grad_plus(j) - grad_minus(j)) / (2.0 * epsilon) - hess(i, j));

                if ((print_level == 1 && fabs(hess_fd_i_j - hess(i, j)) > epsilon) || print_level > 1)
                    std::cout << "Optimization::checkHessian: hess_fd[" << i << ", " << j << "] = " << hess_fd_i_j
                    << ", hess[" << i << ", " << j << "] = " << hess(i, j)
                    << ", error ratio = " << error_ratio << std::endl;

                hess_fd(i, j) = hess_fd_i_j;
            }
        }
        else
        {
            if (print_level > 0)
                std::cout << "Optimization::checkHessian: Gradient evaluation failed!" << std::endl;
            hess_fd.row(i).setConstant(NAN);
        }
    }

    error = hess - hess_fd;
    return error.allFinite() && error.cwiseAbs().maxCoeff() < epsilon;
}
