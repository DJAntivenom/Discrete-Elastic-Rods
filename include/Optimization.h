#ifndef __DER_OPTIMIZATION_H__
#define __DER_OPTIMIZATION_H__

#include <Eigen/Eigen>

template <class data_type = double>
class Optimization
{
public:
    using Float = data_type;
    using SparseMatrix = Eigen::SparseMatrix<Float>;
    using MatrixX = Eigen::MatrixX<Float>;
    using VectorX = Eigen::VectorX<Float>;
    using TripletList = std::vector<Eigen::Triplet<Float>>;
    using Solver = Eigen::SparseLU<SparseMatrix>;

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
    std::function<bool(const VectorX &, Float &)> objective_function;
    std::function<bool(const VectorX &, Float &, VectorX &)> gradient_function;
    std::function<bool(const VectorX &, Float &, VectorX &, TripletList &)> hessian_function;

    Optimization()
    {
        objective_function = [&](const VectorX &y, Float &energy)
        {
            (void)y;
            (void)energy;
            return false;
        };
        gradient_function = [&](const VectorX &y, Float &energy, VectorX &gradient)
        {
            (void)y;
            (void)energy;
            (void)gradient;
            return false;
        };
        hessian_function = [&](const VectorX &y, Float &energy, VectorX &gradient, TripletList &hessian)
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
    bool linearSolve(const TripletList &hessian, const VectorX &b, VectorX &x);

private:
    OptimizationStatus getDirectionGradientDescent(const VectorX &y, VectorX &dy, Float &initial_objective_value);

    OptimizationStatus getDirectionNewton(const VectorX &y, VectorX &dy, Float &initial_objective_value);

public:
    OptimizationStatus step(VectorX &y);

    /// Find step along search direction dy that decreases the objective value. Update y accordingly.
    bool lineSearch(VectorX &y, const VectorX &dy, Float initial_objective_value);

    bool lineSearch(VectorX &y, const VectorX &dy);

    /// Check gradient for scalar function of vector argument.
    static bool checkGradient(const VectorX &y, VectorX &error,
                              const std::function<bool(const VectorX &, Float &)> &func,
                              const std::function<bool(const VectorX &, VectorX &)> &grad_func, Float epsilon,
                              int print_level = 0);

    /// Check hessian for scalar function of vector argument.
    static bool checkHessian(const VectorX &y, MatrixX &error,
                             const std::function<bool(const VectorX &, VectorX &)> &grad_func,
                             const std::function<bool(const VectorX &, MatrixX &)> &hess_func, Float epsilon,
                             int print_level = 0);
};

#include "Optimization.inl"

#endif
