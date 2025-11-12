#ifndef __DER_OPTIMIZATION_H__
#define __DER_OPTIMIZATION_H__

class Optimization
{
public:
    enum class Optimizer
    {
        GRADIENT_DESCENT,
        NEWTON
    };

    enum class OptimizationStatus
    {
        CONVERGED,
        SUCCESS,
        FAILURE
    };

public:
    Optimizer optimizer = Optimizer::GRADIENT_DESCENT;
    int tolerance_exponent = -16;

public:
    /// Gradient and hessian functions also compute objective value and lower-order derivatives.
    std::function<bool(const Eigen::VectorXf &, double &)> objective_function;
    std::function<bool(const Eigen::VectorXf &, double &, Eigen::VectorXf &)> gradient_function;
    std::function<bool(const Eigen::VectorXf &, double &, Eigen::VectorXf &, std::vector<Eigen::Triplet<double>> &)> hessian_function;

    Optimization()
    {
        objective_function = [&](const Eigen::VectorXf &y, double &energy)
        {
            return false;
        };
        gradient_function = [&](const Eigen::VectorXf &y, double &energy, Eigen::VectorXf &gradient)
        {
            return false;
        };
        hessian_function = [&](const Eigen::VectorXf &y, double &energy, Eigen::VectorXf &gradient, std::vector<Eigen::Triplet<double>> &hessian)
        {
            return false;
        };
    }

public:
    /// Solve for x in Ax = b. Return true on success.
    bool linearSolve(const std::vector<Eigen::Triplet<double>> &hessian, const Eigen::VectorXf &b, Eigen::VectorXf &x);

private:
    OptimizationStatus getDirectionGradientDescent(const Eigen::VectorXf &y, Eigen::VectorXf &dy, double &initial_objective_value);

    OptimizationStatus getDirectionNewton(const Eigen::VectorXf &y, Eigen::VectorXf &dy, double &initial_objective_value);

public:
    OptimizationStatus step(Eigen::VectorXf &y);

    /// Find step along search direction dy that decreases the objective value. Update y accordingly.
    bool lineSearch(Eigen::VectorXf &y, const Eigen::VectorXf &dy, double initial_objective_value);

    bool lineSearch(Eigen::VectorXf &y, const Eigen::VectorXf &dy);

    /// Check gradient for scalar function of vector argument.
    static bool checkGradient(const Eigen::VectorXf &y, Eigen::VectorXf &error,
                              const std::function<bool(const Eigen::VectorXf &, double &)> &func,
                              const std::function<bool(const Eigen::VectorXf &, Eigen::VectorXf &)> &grad_func, double epsilon,
                              int print_level = 0);

    /// Check hessian for scalar function of vector argument.
    static bool checkHessian(const Eigen::VectorXf &y, Eigen::MatrixXf &error,
                             const std::function<bool(const Eigen::VectorXf &, Eigen::VectorXf &)> &grad_func,
                             const std::function<bool(const Eigen::VectorXf &, Eigen::MatrixXf &)> &hess_func, double epsilon,
                             int print_level = 0);
};

#endif
