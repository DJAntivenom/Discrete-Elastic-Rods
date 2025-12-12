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
         * We want to compress the rod to 0 < q < 1 of its original length.
         * For simplicity the rod will be bent s.t. there is a single kink.
         * Its new shape combined with the original straight rod compressed
         * to q*length will form an issceles triangle. With pythagoras we
         * can calculate h
         *
         *    q*length
         * ---------------
         * \            /
         *  \          /
         *   \        /
         * l/2\      / l/2
         *     \    /
         *      \  /
         *       \/
         *
         * => h = length * sqrt(1/4 - q²/4)
         *
         * Some care has to be taken if n is even
         *    q*length
         * ---------------
         * \            /
         *  \          /
         *   \        /
         *    \      / nl/(2(n+1))
         *     \    /
         *      \__/
         *      l/(n+1)
         */
        m_rest_omega.setZero();
        const Float q = 0.9;
        m_initial_positions.col(0) = Vector3{ -q * length * 0.5, 0, 0 };
        m_initial_positions.col(n + 1) = Vector3{ q * length * 0.5, 0, 0 };
        const Float l = length / (n + 1);

        Vector3 x0_to_kink, xnp1_to_kink;
        if (n % 2 == 1)
        {
            const Float h = length * std::sqrt(0.25 - q * q * 0.25);
            const Vector3 kink{ 0,-h,0 };

            x0_to_kink = (kink - m_initial_positions.col(0)).normalized();
            xnp1_to_kink = (kink - m_initial_positions.col(n + 1)).normalized();
            m_initial_positions.col((n + 2) / 2) = kink;
        }
        else
        {
            const Float h = std::sqrt(std::pow(n * length / (2 * (n + 1)), 2) - std::pow((length * q - length / (n + 1)) / 2, 2));

            x0_to_kink = (Vector3{ -length * 0.5 / (n + 1), -h, 0 } - m_initial_positions.col(0)).normalized();
            xnp1_to_kink = (Vector3{ length * 0.5 / (n + 1), -h, 0 } - m_initial_positions.col(n + 1)).normalized();
        }

        /**
         * If n is odd we have a vertex at the kink and there are
         * (n+2) / 2 vertices per leg but we can skip the fixed positions.
         * If n is even there is no vertex at the kink but there are still
         * (n+2) / 2 vertices per leg and we can skip the fixed positions
         */
        for (uint64_t i = 1; i < (n + 2) / 2; ++i)
        {
            m_initial_positions.col(i) = m_initial_positions.col(0) + x0_to_kink * l * i;
            m_initial_positions.col(n + 1 - i) = m_initial_positions.col(n + 1) + xnp1_to_kink * l * i;
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
