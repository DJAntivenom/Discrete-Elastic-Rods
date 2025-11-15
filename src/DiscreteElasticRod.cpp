#include <DiscreteElasticRod.hpp>

#include <iostream>

DiscreteElasticRod::DiscreteElasticRod(uint64_t n,
                                       float alpha, float beta,
                                       float radius, float theta_zero, float theta_n) :
    m_vertex_positions(3, (n + 2)),
    m_vertex_velocities(3, (n + 2)),
    m_bishop_frame(3, n + 1),
    m_edge_theta(n + 1),
    m_edge_length(n + 1),
    m_l_i(n + 2),
    m_n(n),
    m_radius(radius),
    m_is_straight_isotropic(true),
    m_alpha(alpha),
    m_beta(beta)
{
    Eigen::Vector3f start{ -0.5f, 0.f, 0.f };
    Eigen::Vector3f increment{ 0.1f, 0.0f, 0.0f };
    for (uint64_t i = 0; i < n + 2; ++i)
    {
        m_vertex_positions.col(i) = start + i * increment;
        if (i > 0)
            m_edge_length[i - 1] = 0.1;
    }
    m_l_i[0] = m_edge_length[0] * 2;
    m_l_i.segment(1, m_n) = m_edge_length.head(m_n) + m_edge_length.tail(m_n);
    m_l_i[m_n + 1] = m_edge_length[m_n] * 2;

    m_total_rod_length = m_edge_length.sum();

    m_vertex_velocities.setZero(); //starting at rest

    m_bishop_frame_vector = Eigen::Vector3f(0.f, 0.f, 0.1f); //perpendicular to first edge
    /**
     * TODO: once we don't just use naturally straigth rods anymore, this has to be done
     * differently.
     */
    m_bishop_frame.colwise() = m_bishop_frame_vector;

    m_edge_theta.setZero();
    m_edge_theta(0) = theta_zero;
    m_edge_theta(n) = theta_n;

    /** TODO: allow anisotropic rods, by minimizing E_bend wrt. B */
    m_B_matrix = Eigen::Vector2f::Constant(alpha).asDiagonal();

    m_w_overbar = getMaterialCurvature(m_edge_theta);
}

void DiscreteElasticRod::update(double delta_time, size_t max_newton_iterations)
{
    /* algorithm outline */

    /// 4., 5. apply torque and integrate rigid body
    /// is handled by just setting the handle to some position

    /// 6., 7. compute forces (given in forumla 11 and above (sec. 7.1))
    /// integrate centerline => apply symplectic euler, see ex.1 handout for formula
    /// initial velocity is zero
    doSymplecticEuler(delta_time);

    /// 8. TODO: enforce constraints to guarantee inextensibility

    /// 9. TODO: handle collisions (not discussed in detail in paper, but reference in paper)

    /// 10. update natural bishop frame, i.e. apply rotation (P_i in paper)
    transportBishopFrame();

    /// 11. quasistatic material frame update (Newton according to equation 4 in paper)
    /// => Use newton solver from exercises
    applyTwist(max_newton_iterations);
}

void DiscreteElasticRod::setVertexPosition(uint64_t vertex_index, const Eigen::Vector3f &new_position)
{
    m_vertex_positions.col(vertex_index) = new_position;

    m_is_straight_isotropic = false;
    transportBishopFrame();
}

void DiscreteElasticRod::randomizeVertexPositions()
{
    m_vertex_positions.block(1, 0, 2, m_n + 2).setRandom() *= 0.1f;

    m_is_straight_isotropic = false;
    transportBishopFrame();
}

static Eigen::Vector3f hsvToRgb(const Eigen::Vector3f &hsv)
{
    const float hue = hsv.x();
    const float sat = hsv.y();
    const float val = hsv.z();

    assert(sat >= 0.f && sat <= 1.f);
    assert(val >= 0.f && val <= 1.f);

    /* convert hue to degree and get first equivalent value */
    const float hue_sextant = fmod(hue * 180.f * M_1_PI, 360.f) / 60.f;
    const int64_t hue_integer = static_cast<int64_t>(hue_sextant);
    const float hue_fraction = hue_sextant - hue_integer;

    const float p = val * (1.0 - sat);
    const float q = val * (1.0 - (sat * hue_fraction));
    const float t = val * (1.0 - (sat * (1.0 - hue_fraction)));

    switch (hue_integer)
    {
    case 0:
        return Eigen::Vector3f{ val, t, p };
    case 1:
        return Eigen::Vector3f{ q, val, p };
    case 2:
        return Eigen::Vector3f{ p, val, t };
    case 3:
        return Eigen::Vector3f{ p, q, val };
    case 4:
        return Eigen::Vector3f{ t, p, val };
    case 5:
    default:
        return Eigen::Vector3f{ val, p, q };
    }
}

polyscope::SurfaceMesh *DiscreteElasticRod::registerSurfaceMesh(const std::string &name,
                                                                uint32_t vertices_per_ring) const
{
    using namespace Eigen;
    constexpr float saturation = 1.f;
    constexpr float value = 1.f;

    /**
     * We have the starting vertex per edge, there are
     * (n + 1) edges, there is a single anchor at the end of the rod
     */
    const uint32_t number_of_centerline_anchors = (m_n + 1) + 1;

    /**
     * In addition to the `vertices_per_ring` we need for the rings surrounding the
     * anchors, we also need two vertices at the end and beginning of the rod
     */
    const uint32_t number_of_vertices = number_of_centerline_anchors * vertices_per_ring + 2;

    /**
     * Each ring-segment is connected to the next ring with a rectangle, which is 2 triangles.
     * Additionally, the ends of the rod are connected to the ring.
     */
    const uint32_t number_of_faces = (number_of_centerline_anchors - 1) * vertices_per_ring * 2 + 2 * vertices_per_ring;

    MatrixX3f colors(number_of_vertices, 3);
    MatrixX3f vertices(number_of_vertices, 3);
    MatrixX3i faces(number_of_faces, 3);

    const auto tangents = getTangents();

    const float inv_vertices_per_ring = 1.f / vertices_per_ring;

    /**
     * We go through all edges and set ring vertices
     */
    uint32_t vertex_index = 0, face_index = 0;
    for (uint32_t edge_index = 0; edge_index <= m_n; ++edge_index)
    {
        const Vector3f anchor_position = m_vertex_positions.col(edge_index);
        /* axis of rotation is length-weighted average of vertex-adjacent edges */
        const Vector3f axis_of_rotation = (edge_index > 0
                                           ? (tangents.col(edge_index) * m_edge_length[edge_index]
                                              + tangents.col(edge_index - 1) * m_edge_length[edge_index - 1]).normalized()
                                           : tangents.col(edge_index));

        for (uint32_t ring_index = 0; ring_index < vertices_per_ring; ++ring_index, ++vertex_index)
        {
            /* vertex positions are rotated about anchors, */
            const float ring_angle = (2 * M_PI) * ring_index * inv_vertices_per_ring;
            const float theta = (edge_index > 0
                                 ? (m_edge_theta[edge_index] * m_edge_length[edge_index]
                                    + m_edge_theta[edge_index - 1] * m_edge_length[edge_index - 1])
                                 / (m_edge_length[edge_index] + m_edge_length[edge_index - 1])
                                 : m_edge_theta[edge_index]);

            /* handles rotate clockwise */
            const float angle = ring_angle - theta;

            vertices.row(vertex_index).transpose() = anchor_position +
                AngleAxisf(angle, axis_of_rotation) * m_bishop_frame.col(edge_index) * m_radius;
            colors.row(vertex_index) = hsvToRgb({ angle, saturation, value });

            /* faces are connected to previous ring */
            if (edge_index > 0)
            {
                const uint32_t next_vertex = (ring_index == vertices_per_ring - 1
                                              ? vertex_index + 1 - vertices_per_ring
                                              : vertex_index + 1);
                faces.row(face_index++) <<
                    vertex_index,
                    next_vertex,
                    vertex_index - vertices_per_ring;
                faces.row(face_index++) <<
                    next_vertex,
                    next_vertex - vertices_per_ring,
                    vertex_index - vertices_per_ring;
            }
        }
    }

    /* do last anchor */
    for (uint32_t ring_index = 0; ring_index < vertices_per_ring; ++ring_index, ++vertex_index)
    {
        /* vertex positions are rotated about anchors */
        const float angle = (2 * M_PI) * ring_index * inv_vertices_per_ring
            - m_edge_theta[m_n];
        vertices.row(vertex_index).transpose() = m_vertex_positions.col(m_n + 1) +
            AngleAxisf(angle, tangents.col(m_n)) * m_bishop_frame.col(m_n) * m_radius;
        colors.row(vertex_index) = hsvToRgb({ angle, saturation, value });

        /* faces are connected to previous ring */
        const uint32_t next_vertex = (ring_index == vertices_per_ring - 1
                                      ? vertex_index + 1 - vertices_per_ring
                                      : vertex_index + 1);
        faces.row(face_index++) <<
            vertex_index,
            next_vertex,
            vertex_index - vertices_per_ring;
        faces.row(face_index++) <<
            next_vertex,
            next_vertex - vertices_per_ring,
            vertex_index - vertices_per_ring;
    }

    /* set end of rods */
    vertices.row(vertex_index++) = m_vertex_positions.col(0).transpose();
    vertices.row(vertex_index++) = m_vertex_positions.col(m_n + 1).transpose();
    colors.bottomRows(2).setOnes();
    assert(vertex_index == number_of_vertices);

    /* do ends of rod faces */
    for (uint32_t ring_index = 0; ring_index < vertices_per_ring; ++ring_index)
    {
        faces.row(face_index++) <<
            ring_index,
            (ring_index + 1) % vertices_per_ring,
            number_of_vertices - 2;

        const uint32_t vertex_index =
            number_of_vertices - 2 - vertices_per_ring + ring_index;
        const uint32_t next_vertex = (ring_index == vertices_per_ring - 1
                                      ? vertex_index + 1 - vertices_per_ring
                                      : vertex_index + 1);
        faces.row(face_index++) <<
            vertex_index,
            next_vertex,
            number_of_vertices - 1;
    }
    assert(face_index == number_of_faces);

    auto mesh = polyscope::registerSurfaceMesh(name, vertices, faces);

    mesh->addVertexColorQuantity("rotation color", colors)->setEnabled(true);

    return mesh;
}

Eigen::Matrix4Xf DiscreteElasticRod::getMaterialCurvature(const Eigen::VectorXf &theta) const
{
    using namespace Eigen;

    Matrix4Xf curvature(4, m_n + 2);
    const auto binormals = getBinormals();
    const auto tangents = getTangents();

    for (uint32_t vertex = 1; vertex <= m_n; ++vertex)
    {
        /** NOTE: need to subtract 1 in binormals.col because binormals is just inner vertices */
        {
            const uint32_t j = vertex - 1;
            const Vector3f m_1_i_min_1 = AngleAxisf(theta[j], tangents.col(j)) * m_bishop_frame.col(j);
            const Vector3f m_2_i_min_1 = AngleAxisf(theta[j] + M_PI_2, tangents.col(j)) * m_bishop_frame.col(j);
            curvature.block<2, 1>(0, vertex) <<
                binormals.col(vertex - 1).dot(m_2_i_min_1),
                binormals.col(vertex - 1).dot(m_1_i_min_1);
        }

        {
            const uint32_t j = vertex;
            const Vector3f m_1_i = AngleAxisf(theta[j], tangents.col(j)) * m_bishop_frame.col(j);
            const Vector3f m_2_i = AngleAxisf(theta[j] + M_PI_2, tangents.col(j)) * m_bishop_frame.col(j);
            curvature.block<2, 1>(2, vertex) <<
                binormals.col(vertex - 1).dot(m_2_i),
                binormals.col(vertex - 1).dot(m_1_i);
        }
    }

    return curvature;
}
