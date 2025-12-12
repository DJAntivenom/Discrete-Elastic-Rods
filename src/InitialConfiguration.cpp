#include <InitialConfiguration.hpp>

class IsotropicStraightRest : public InitialConfiguration
{
public:
    IsotropicStraightRest(uint64_t n, float radius, Float length) :
        InitialConfiguration(n, radius)
    {
        //const Float length = 1.0f;
        const Vector3 start{ -0.5f * length, 0, 0 };
        const Vector3 increment{ length / (n + 1), 0, 0 };

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
        /**
         * We want to have the rod form a circle-segment defined by angle q using
         * equally spaced angle segments.
         */
        m_rest_omega.setZero();
        /* 30° */
        constexpr double q = M_PI / 6;
        const double h = length / (2 * std::tan(q * 0.5));
        const double r = std::sqrt(h * h + length * length * 0.25);
        const Vector3 center{ 0, h, 0 };

        for (uint64_t i = 0; i <= n + 1; ++i)
        {
            const Float angle = q * i / (n + 1) + (M_PI + M_PI_2 - q / 2);
            m_initial_positions.col(i) = r * Vector3{ std::cos(angle), std::sin(angle), 0 } + center;
        }
    }
};

std::unique_ptr<InitialConfiguration> getInitialConfiguration(InitialConfiguration::ConfigType type,
                                                              uint64_t n,
                                                              float radius,
                                                              Float length)
{
    assert(length > 0);

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
