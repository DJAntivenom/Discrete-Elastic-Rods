#include <DiscreteElasticRod.hpp>

#include <iostream>

#include "common.h"

DiscreteElasticRod::DiscreteElasticRod(const InitialConfiguration &ic,
                                       Float alpha, Float beta) :
    m_vertex_positions(ic.getInitialPositions()),
    m_vertex_velocities(3, ic.getN() + 2),
    m_bishop_frame_vector{ 0., 0., 1. },
    m_bishop_frame(3, ic.getN() + 1),
    m_edge_theta(ic.getN() + 1),
    m_edge_length(ic.getN() + 1),
    m_vertex_mass(ic.getN() + 2),
    m_l_i(ic.getN() + 2),
    m_n(ic.getN()),
    m_radius(ic.getRadius()),
    m_is_straight_isotropic(ic.isStraight() && ic.isIsotropic()),
    m_alpha(alpha),
    m_beta(beta)
{
    for (uint64_t i = 1; i < m_n + 2; ++i)
    {
        m_edge_length[i - 1] = (m_vertex_positions.col(i) - m_vertex_positions.col(i - 1)).norm();
    }
    m_l_i[0] = m_edge_length[0];
    m_l_i.segment(1, m_n) = m_edge_length.head(m_n) + m_edge_length.tail(m_n);
    m_l_i[m_n + 1] = m_edge_length[m_n];

    m_vertex_mass.setConstant(1.0f);

    m_total_rod_length = m_edge_length.sum();

    m_vertex_velocities.setZero(); //starting at rest

    m_edge_theta.setZero();

    /** TODO: allow anisotropic rods, by minimizing E_bend wrt. B */
    m_B_matrix = Vector2::Constant(alpha).asDiagonal();

    m_w_overbar = ic.getRestOmega();

    transportBishopFrame();
}

/*
 * !ONLY TO BE USED FOR THE CUTTING METHOD DISCRETEELASTICROD::CUTATVERTEX!
 */
DiscreteElasticRod::DiscreteElasticRod(Matrix3X &vertex_positions,
                                        VectorX &edge_lengths,
                                        Matrix3X &vertex_velocities,
                                        Vector3 &bishop_frame_vector,
                                        Matrix3X &bishop_frame,
                                        VectorX &edge_theta,
                                        VectorX &vertex_mass,
                                        VectorX &l_i,
                                        Matrix2 &B_matrix,
                                        Matrix4X &w_overbar,
                                        int n,
                                        Float radius,
                                        bool is_straight_isotropic,
                                        Float alpha,
                                        Float beta) :
    m_vertex_positions(vertex_positions),
    m_edge_length(edge_lengths),
    m_vertex_velocities(vertex_velocities),
    m_bishop_frame_vector(bishop_frame_vector),
    m_bishop_frame(bishop_frame),
    m_edge_theta(edge_theta),
    m_vertex_mass(vertex_mass),
    m_l_i(l_i),
    m_n(n),
    m_radius(radius),
    m_is_straight_isotropic(is_straight_isotropic),
    m_alpha(alpha),
    m_beta(beta),
    m_B_matrix(B_matrix),
    m_w_overbar(w_overbar)
{
    m_total_rod_length = m_edge_length.sum();
}

void DiscreteElasticRod::update(double delta_time, size_t max_newton_iterations)
{
#ifdef KEEP_TURNING
    m_edge_theta[m_n] += 10 * M_PI * delta_time;
#endif
    /* algorithm outline */

    /// 4., 5. apply torque and integrate rigid body
    /// is handled by just setting the handle to some position

    /// 6., 7. compute forces (given in forumla 11 and above (sec. 7.1))
    /// integrate centerline => apply symplectic euler, see ex.1 handout for formula
    /// initial velocity is zero
    doSymplecticEuler(delta_time);

    /// 8. TODO: enforce constraints to guarantee inextensibility
    applyConstraints(max_newton_iterations);
    print_debug("post constraints");
    /// 9. TODO: handle collisions (not discussed in detail in paper, but reference in paper)

    /// 10. update natural bishop frame, i.e. apply rotation (P_i in paper)
    transportBishopFrame();
    print_debug("post transporting bishop frame");
    print_debug(m_vertex_positions);

    /// 11. quasistatic material frame update (Newton according to equation 4 in paper)
    /// => Use newton solver from exercises
    applyTwist(max_newton_iterations);
    print_debug("post-twist");
    print_debug(m_vertex_positions);
}

void DiscreteElasticRod::setVertexPosition(uint64_t vertex_index, const Vector3 &new_position)
{
    m_vertex_positions.col(vertex_index) = new_position;

    transportBishopFrame();
}

void DiscreteElasticRod::randomizeVertexPositions()
{
    m_vertex_positions.block(1, 0, 2, m_n + 2).setRandom() *= 0.1f;

    m_edge_length = getEdges().colwise().norm();

    m_l_i[0] = m_edge_length[0];
    m_l_i.segment(1, m_n) = m_edge_length.head(m_n) + m_edge_length.tail(m_n);
    m_l_i[m_n + 1] = m_edge_length[m_n];

    m_total_rod_length = m_edge_length.sum();

    m_vertex_velocities.setZero(); //starting at rest

    m_is_straight_isotropic = false;
    transportBishopFrame();
}

template <class Vector3>
static Vector3 hsvToRgb(const Vector3 &hsv)
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
        return Vector3{ val, t, p };
    case 1:
        return Vector3{ q, val, p };
    case 2:
        return Vector3{ p, val, t };
    case 3:
        return Vector3{ p, q, val };
    case 4:
        return Vector3{ t, p, val };
    case 5:
    default:
        return Vector3{ val, p, q };
    }
}

polyscope::SurfaceMesh *DiscreteElasticRod::registerSurfaceMesh(const std::string &name,
                                                                uint32_t vertices_per_ring) const
{
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

    MatrixX3 colors(number_of_vertices, 3);
    MatrixX3 vertices(number_of_vertices, 3);
    Eigen::MatrixX3i faces(number_of_faces, 3);

    const auto tangents = getTangents();

    const float inv_vertices_per_ring = 1.f / vertices_per_ring;

    /**
     * We go through all edges and set ring vertices
     */
    uint32_t vertex_index = 0, face_index = 0;
    for (uint32_t edge_index = 0; edge_index <= m_n; ++edge_index)
    {
        const Vector3 anchor_position = m_vertex_positions.col(edge_index);
        /* axis of rotation is length-weighted average of vertex-adjacent edges */
        const Vector3 axis_of_rotation = (edge_index > 0
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
                AngleAxis(angle, axis_of_rotation) * m_bishop_frame.col(edge_index) * m_radius;
            colors.row(vertex_index) = hsvToRgb<Vector3>({ angle, saturation, value });

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
            AngleAxis(angle, tangents.col(m_n)) * m_bishop_frame.col(m_n) * m_radius;
        colors.row(vertex_index) = hsvToRgb<Vector3>({ angle, saturation, value });

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

Matrix4X DiscreteElasticRod::getMaterialCurvature(const VectorX &theta) const
{
    Matrix4X curvature(4, m_n + 2);
    const auto binormals = getBinormals();
    const auto tangents = getTangents();

    for (uint32_t vertex = 1; vertex <= m_n; ++vertex)
    {
        /** NOTE: need to subtract 1 in binormals.col because binormals is just inner vertices */
        {
            const uint32_t j = vertex - 1;
            const Vector3 m_1_i_min_1 = AngleAxis(theta[j], tangents.col(j)) * m_bishop_frame.col(j);
            const Vector3 m_2_i_min_1 = AngleAxis(theta[j] + M_PI_2, tangents.col(j)) * m_bishop_frame.col(j);
            curvature.block<2, 1>(0, vertex) <<
                binormals.col(vertex - 1).dot(m_2_i_min_1),
                binormals.col(vertex - 1).dot(m_1_i_min_1);
        }

        {
            const uint32_t j = vertex;
            const Vector3 m_1_i = AngleAxis(theta[j], tangents.col(j)) * m_bishop_frame.col(j);
            const Vector3 m_2_i = AngleAxis(theta[j] + M_PI_2, tangents.col(j)) * m_bishop_frame.col(j);
            curvature.block<2, 1>(2, vertex) <<
                binormals.col(vertex - 1).dot(m_2_i),
                binormals.col(vertex - 1).dot(m_1_i);
        }
    }

    return curvature;
}

std::pair<DiscreteElasticRod, DiscreteElasticRod> DiscreteElasticRod::cutAtVertex(int i) {
    assert(i > 0 && i < m_n + 1);
    print_debug(m_n + 2);
    print_debug(i);

    Matrix3X vertex_positions_1 = m_vertex_positions.leftCols(i + 1);
    Matrix3X vertex_positions_2 = m_vertex_positions.rightCols(m_n + 2 - i);

    print_debug(m_vertex_positions.cols());
    print_debug(vertex_positions_1.cols());
    print_debug(vertex_positions_2.cols());

    VectorX edge_lengths_1 = m_edge_length.head(i);
    VectorX edge_lengths_2 = m_edge_length.tail(m_n + 1 - i);

    print_debug(m_edge_length.size());
    print_debug(edge_lengths_1.size());
    print_debug(edge_lengths_2.size());

    Matrix3X vertex_velocities_1 = m_vertex_velocities.leftCols(i + 1);
    Matrix3X vertex_velocities_2 = m_vertex_velocities.rightCols(m_n + 2 - i);

    print_debug(m_vertex_velocities.cols());
    print_debug(vertex_velocities_1.cols());
    print_debug(vertex_velocities_2.cols());

    Vector3 bishop_frame_vector_1 = m_bishop_frame_vector;
    Vector3 bishop_frame_vector_2 = m_bishop_frame_vector;

    Matrix3X bishop_frame_1 = m_bishop_frame.leftCols(i);
    Matrix3X bishop_frame_2 = m_bishop_frame.rightCols(m_n + 1 - i);

    VectorX edge_theta_1 = m_edge_theta.head(i);
    VectorX edge_theta_2 = m_edge_theta.tail(m_n + 1 - i);

    VectorX vertex_mass_1 = m_vertex_mass.head(i + 1);
    VectorX vertex_mass_2 = m_vertex_mass.tail(m_n + 2 - i);

    VectorX l_i_1 = m_l_i.head(i + 1);
    VectorX l_i_2 = m_l_i.tail(m_n + 2 - i);

    Matrix4X w_overbar_1 = m_w_overbar.leftCols(i + 1);
    Matrix4X w_overbar_2 = m_w_overbar.rightCols(m_n + 2 - i);

    int n_1 = i - 1;
    int n_2 = m_n - i;

    print_debug(m_n);
    print_debug(n_1);
    print_debug(n_2);

    DiscreteElasticRod discrete_elastic_rod_1(vertex_positions_1,
        edge_lengths_1, vertex_positions_1, bishop_frame_vector_1,
        bishop_frame_1, edge_theta_1, vertex_mass_1,
        l_i_1, m_B_matrix, w_overbar_1, n_1,
        m_radius, m_is_straight_isotropic, m_alpha, m_beta) ;

    DiscreteElasticRod discrete_elastic_rod_2(vertex_positions_2,
        edge_lengths_2, vertex_positions_2, bishop_frame_vector_2,
        bishop_frame_2, edge_theta_2, vertex_mass_2,
        l_i_2, m_B_matrix, w_overbar_2, n_2,
        m_radius, m_is_straight_isotropic, m_alpha, m_beta);

    //TODO: Test that no data is duplicated (except for the final vertex) or lost this way
    //TODO: how to delete this rod after creating the other two, and make sure that these new ones are animated?
    return std::make_pair(discrete_elastic_rod_1, discrete_elastic_rod_2);
}
