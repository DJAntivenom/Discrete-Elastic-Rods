#include <InitialConfiguration.hpp>

class IsotropicStraightRest : public InitialConfiguration
{
public:
    IsotropicStraightRest(uint64_t n, float radius, Float length) :
        InitialConfiguration(n, radius)
    {
        //const Float length = 1.0f;
        const Vector3 start{ -0.5f * length, 0, 0 };
        const Vector3 increment{ length / (n + 2), 0, 0 };

        for (uint32_t i = 0; i < n + 2; ++i)
        {
            m_initial_positions.col(i) = start + i * increment;
        }

        m_rest_omega.setZero();
    }
};

class IsotropicStraightUnderPressure : public InitialConfiguration
{
public:
    IsotropicStraightUnderPressure(const uint64_t n, const float radius, const Float length) :
        InitialConfiguration(n, radius)
    {
        const Float horizontal_ellipse_radius = length / 2;
        const Float vertical_ellipse_radius = length / 20;

        for (uint32_t i = 0; i < n + 2; ++i)
        {
            const Float t = (1 - static_cast<Float>(i) / (n + 1)) * M_PI;
            m_initial_positions.col(i) = Vector3{ horizontal_ellipse_radius * std::cos(t), -vertical_ellipse_radius * std::sin(t), 0 };
        }

        m_rest_omega.setZero();
    }
};

std::unique_ptr<InitialConfiguration> getInitialConfiguration(InitialConfiguration::ConfigType type,
                                                              uint64_t n,
                                                              float radius,
                                                              Float length)
{
    assert(length > 0 && length <= (n + 1) / 10);

    switch (type)
    {
    case InitialConfiguration::STRAIGHT_ISOTROPIC_AT_REST:
        return std::make_unique<IsotropicStraightRest>(n, radius, length);
    case InitialConfiguration::STRAIGHT_ISOTROPIC_PRESSURE:
        return std::make_unique<IsotropicStraightUnderPressure>(n, radius, length);
    case InitialConfiguration::UNKNOWN:
        throw std::runtime_error("Unknown configuration type");
    default:
        throw std::runtime_error("Initial Configuration type not implemented!");
    }
}
