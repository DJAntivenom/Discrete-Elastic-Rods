#include <DiscreteElasticRod.hpp>

#include <stdexcept>

void DiscreteElasticRod::applyTwist(size_t max_newton_iterations)
{
    (void)max_newton_iterations;
    throw std::runtime_error(std::string("DiscreteElasticRod::applyTwist(): Not implemented"));
}
