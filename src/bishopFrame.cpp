#include <DiscreteElasticRod.hpp>

#include <algorithm>
#include <stdexcept>

void DiscreteElasticRod::transportBishopFrame()
{
    const auto tangents = getTangents();
    const auto binormals = getBinormals();

    /* reproject initial u_0 to be perpendicular to t_0 */
    m_bishop_frame_vector -= tangents.col(0) * tangents.col(0).dot(m_bishop_frame_vector);
    m_bishop_frame_vector.normalize();

    m_bishop_frame.col(0) = m_bishop_frame_vector;

    for (uint32_t i = 0; i < m_n; ++i)
    {
        if (binormals.col(i).isZero())
        {
            m_bishop_frame.col(i + 1) = m_bishop_frame.col(i);
        }
        else
        {
            const Float angleCos = tangents.col(i).dot(tangents.col(i + 1));
            const Float angle = std::acos(std::clamp(angleCos, -1.0, 1.0));
            IsometricTransform3D rotation(AngleAxis(angle, binormals.col(i).normalized()));

            m_bishop_frame.col(i + 1) = rotation.linear() * m_bishop_frame.col(i);
        }
    }
}

DiscreteElasticRod::Matrix3 DiscreteElasticRod::getInterpolatedBishopFrame(double alpha)
{
    const float interpolated_length = m_total_rod_length * alpha;

    uint64_t edge_index = 0;
    float rest_length = interpolated_length;
    for (; rest_length > m_edge_length[edge_index]; ++edge_index)
        rest_length -= m_edge_length[edge_index];

    const auto tangents = getTangents();
    Matrix3 result;
    result.row(0) = (m_vertex_positions.col(edge_index) + tangents.col(edge_index) * rest_length).transpose();
    result.row(1) = m_bishop_frame.col(edge_index).transpose();
    result.row(2) = m_bishop_frame.col(edge_index).cross(tangents.col(edge_index)).transpose();

    return result;
}
