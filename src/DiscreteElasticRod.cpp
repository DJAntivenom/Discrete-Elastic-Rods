#include <DiscreteElasticRod.hpp>

DiscreteElasticRod::DiscreteElasticRod(uint64_t n) : m_vertex_positions(3 * (n + 2)), m_num_vertices(n + 2)
{
    Eigen::Vector3f start{0.f, 0.5f, 0.f};
    Eigen::Vector3f increment{0.1f, 0.0f, 0.0f};
    for (uint64_t i = 0; i < n + 2; ++i)
    {
        m_vertex_positions.segment<3>(i * 3) = start + i * increment;
    }
}

void DiscreteElasticRod::update(double delta_time)
{
    static float base_time = 0;
    const float y = 0.5f + std::sin(base_time) * 0.1f;
    base_time += delta_time;

    getVertex(m_num_vertices / 2).y() = y;
}
