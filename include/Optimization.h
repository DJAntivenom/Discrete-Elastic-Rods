#ifndef __DER_OPTIMIZATION_H__
#define __DER_OPTIMIZATION_H__

#include <Eigen/Eigen>

class Optimization
{
public:
    using SparseMatrixF = Eigen::SparseMatrix<float>;
    using MatrixXf = Eigen::MatrixXf;
    using VectorXf = Eigen::VectorXf;
    using TripletListF = std::vector<Eigen::Triplet<float>>;
    using Solver = Eigen::SparseLU<SparseMatrixF>;

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
    std::function<bool(const VectorXf &, double &)> objective_function;
    std::function<bool(const VectorXf &, double &, VectorXf &)> gradient_function;
    std::function<bool(const VectorXf &, double &, VectorXf &, TripletListF &)> hessian_function;

    Optimization()
    {
        objective_function = [&](const VectorXf &y, double &energy)
        {
            (void)y;
            (void)energy;
            return false;
        };
        gradient_function = [&](const VectorXf &y, double &energy, VectorXf &gradient)
        {
            (void)y;
            (void)energy;
            (void)gradient;
            return false;
        };
        hessian_function = [&](const VectorXf &y, double &energy, VectorXf &gradient, TripletListF &hessian)
        {
            (void)y;
            (void)energy;
            (void)gradient;
            (void)hessian;
            return false;
        };
    }

public:
    /// Solve for x in Ax = b. Return true on success.
    bool linearSolve(const TripletListF &hessian, const VectorXf &b, VectorXf &x);

private:
    OptimizationStatus getDirectionGradientDescent(const VectorXf &y, VectorXf &dy, double &initial_objective_value);

    OptimizationStatus getDirectionNewton(const VectorXf &y, VectorXf &dy, double &initial_objective_value);

public:
    OptimizationStatus step(VectorXf &y);

    /// Find step along search direction dy that decreases the objective value. Update y accordingly.
    bool lineSearch(VectorXf &y, const VectorXf &dy, double initial_objective_value);

    bool lineSearch(VectorXf &y, const VectorXf &dy);

    /// Check gradient for scalar function of vector argument.
    static bool checkGradient(const VectorXf &y, VectorXf &error,
                              const std::function<bool(const VectorXf &, double &)> &func,
                              const std::function<bool(const VectorXf &, VectorXf &)> &grad_func, double epsilon,
                              int print_level = 0);

    /// Check hessian for scalar function of vector argument.
    static bool checkHessian(const VectorXf &y, MatrixXf &error,
                             const std::function<bool(const VectorXf &, VectorXf &)> &grad_func,
                             const std::function<bool(const VectorXf &, MatrixXf &)> &hess_func, double epsilon,
                             int print_level = 0);
};

#endif
