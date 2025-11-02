#include <DiscreteElasticRod.hpp>

DiscreteElasticRod::DiscreteElasticRod(uint64_t n, bool clamped, float theta_zero, float theta_n) : 
    m_vertex_positions(3 * (n + 2)), m_vertex_velocities(3 * (n + 2)), m_edge_theta(n), m_num_vertices(n + 2), is_clamped(clamped)
{
    Eigen::Vector3f start{0.f, 0.5f, 0.f};
    Eigen::Vector3f increment{0.1f, 0.0f, 0.0f};
    for (uint64_t i = 0; i < n + 2; ++i)
    {
        m_vertex_positions.segment<3>(i * 3) = start + i * increment;
    }
    bishop_frame_vector = Eigen::Vector3f(0.f, 0.f, 0.1f); //perpendicular to first edge
    
    m_vertex_velocities = Eigen::VectorXf::Zero(3 * (n + 2)); //starting at rest

    m_edge_theta = Eigen::VectorXf::Zero(n); 
    if (is_clamped) {
        m_edge_theta(0) = theta_zero;
        m_edge_theta(n) = theta_n;
    }
}

void DiscreteElasticRod::update(double delta_time)
{
    static float base_time = 0;
    const float y = 0.5f + std::sin(base_time) * 0.1f;
    base_time += delta_time;

    getVertex(m_num_vertices / 2).y() = y;
}
