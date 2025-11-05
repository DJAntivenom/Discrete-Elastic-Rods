#include <DiscreteElasticRod.hpp>

#include <stdexcept>

void DiscreteElasticRod::doSymplecticEuler(double delta_time)
{
    (void)delta_time;
    throw std::runtime_error(std::string("DiscreteElasticRod::doSymplecticEuler(double): Not implemented"));
}
