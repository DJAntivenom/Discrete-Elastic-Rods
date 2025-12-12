#ifndef __DER_INITIAL_CONFIGURATION_HPP__
#define __DER_INITIAL_CONFIGURATION_HPP__

#include "common.h"

#include <cstdint>
#include <memory>

class InitialConfiguration
{
public:
    /**
     * \brief Enumeration type fot the different existing initial configurations.
     */
    enum ConfigType
    {
        STRAIGHT_ISOTROPIC_AT_REST,  //< A straight isotropic rod that is unperturbed.
        STRAIGHT_ISOTROPIC_PRESSURE, //< A straight isotropic rod that is pushed together to form a half-circle.
        STRAIGHT_ISOTROPIC_TREFOIL,
        CURVED_ISOTROPIC,
        UNKNOWN,
        CONFIG_COUNT = UNKNOWN
    };

    /**
     * @brief Namey of the existing initial configurations.
     */
    static constexpr const char *const CONFIG_NAMES[] = {
        "Straight Isotropic unperturbed",
        "Straight Isotropic under pressure",
        "Straight Isotropic in trefoil",
        "Curved Isotropic"
    };
    static_assert(sizeof(CONFIG_NAMES) / sizeof(CONFIG_NAMES[0]) == ConfigType::CONFIG_COUNT);

    /// Objects should only be created through the factory function below
    InitialConfiguration() = delete;
    InitialConfiguration(const InitialConfiguration &) = delete;
    InitialConfiguration(InitialConfiguration &&) = delete;
    InitialConfiguration &operator=(const InitialConfiguration &) = delete;
    InitialConfiguration &operator=(InitialConfiguration &&) = delete;
    virtual ~InitialConfiguration() = default;

    /**
     * @brief Whether the rod is straight when in rest.
     * @return `true` iff the rod is straight in rest.
     */
    virtual bool isStraight() const
    {
        return true;
    }

    /**
     * @brief Whether the rod is straight when in rest.
     * @return `true` iff the rod is straight in rest.
     */
    virtual bool isIsotropic() const
    {
        return true;
    }

    /**
     * @brief Get the number of inner vertices.
     * @return The returned value does not count the boundary vertices.
     */
    uint64_t getN() const
    {
        return m_n;
    }

    /**
     * @brief Get the radius of the rod.
     */
    float getRadius() const
    {
        return m_radius;
    }

    /**
     * @brief The initial positions might not be in rest, however, their edge's lengths are the rest length.
     */
    virtual Matrix3X getInitialPositions() const
    {
        return m_initial_positions;
    }

    /**
     * @brief The material curvature at rest. This is omega bar in the paper.
     */
    virtual Matrix4X getRestOmega() const
    {
        return m_rest_omega;
    }

protected:
    InitialConfiguration(uint64_t n, float radius) :
        m_n(n),
        m_radius(radius),
        m_initial_positions(3, n + 2),
        m_rest_omega(4, n + 2)
    {
    };

    const uint64_t m_n;
    const float m_radius;

    /**
     * @brief dim 3x(n+2)
     */
    Matrix3X m_initial_positions;
    /**
     * \brief 4 x (n + 2), every column is
     * {omega_i^[i-1].x, omega_i^[i-1].y, omega_i^[i].x, omega_i^[i].y}.
     */
    Matrix4X m_rest_omega;
};

/**
 * @brief Creates an initial configuration based on the given type.
 * @param type Type of the rod to create. \see InitialConfiguration::ConfigType
 * @param n The number of inner vertices.
 * @param radius The radius of the rod.
 * @param length_between_boundaries Used by some implementations as the initial distance
 * of the boundary-vertices. Must be greater than 0.
 * @return A unique pointer to the created initial configuration.
 */
std::unique_ptr<InitialConfiguration> getInitialConfiguration(InitialConfiguration::ConfigType type,
                                                              uint64_t n,
                                                              float radius,
                                                              Float length_between_boundaries = 1.);

#endif
