#include <DiscreteElasticRod.hpp>
#include <common.h>

static VectorX minimum_distance(const int i, const int j, const Matrix3X &vertices_i, const Matrix3X &vertices_j) {
    Vector3 xi = vertices_i.col(i);
    Vector3 xi_plus = vertices_i.col(i+1);
    Vector3 xj = vertices_j.col(j);
    Vector3 xj_plus = vertices_j.col(j+1);

    Vector3 ei = xi_plus - xi;
    Vector3 ej = xj_plus - xj;

    Vector3 n = ei.cross(ej);
    if(n.squaredNorm() < 1e-8) {
        return VectorX::Constant(5, -1);
    }

    Float t1 = std::clamp<Float>((ej.cross(n).dot(xj - xi)) / n.squaredNorm(), 0, 1);
    Float t2 = std::clamp<Float>((ei.cross(n).dot(xj - xi)) / n.squaredNorm(), 0, 1);

    Vector3 p1 = xi + t1 * ei;
    Vector3 p2 = xj + t2 * ej;

    VectorX res(5);
    res.setZero();
    res.head(3) = (p2 - p1);
    res[3] = t1;
    res[4] = t2;

    return res;
}

void DiscreteElasticRod::calculateContactForces(const int max_iter, const Float d, DiscreteElasticRod &other, const bool same, const Float delta_time) {
    Float err_tol = 1e-8;

    VectorX err_vec = VectorX::Constant((m_vertex_positions.cols() - 1) * (other.m_vertex_velocities.cols() - 1), 1.f);
    for (int k = 0; k < max_iter && err_vec.cwiseGreaterOrEqual(err_tol).any(); ++k) {
        for (int i = 0; i < m_vertex_positions.cols() - 1; i++) {
            for (int j = 0; j < other.m_vertex_positions.cols() - 1; j++) {
                if(err_vec[i * (other.m_vertex_positions.cols() - 1) + j] > 0 && (!same || abs(i - j) > 2)) {
                    auto res = minimum_distance(i, j, m_vertex_positions, other.m_vertex_positions);
                    Vector3 n_ij = res.head(3);
                    Float md_ij = n_ij.norm();
                    if (res[3] != -1 && md_ij < d) {
                        print_debug(res[3]);
                        print_debug("md_ij = " + std::to_string(md_ij));
                        //VectorX res = minimum_distance(i, j, m_vertex_positions, other.m_vertex_positions);
                        print_debug(std::to_string(i) + " and " + std::to_string(j) + " are too close to each other because the distance is " + std::to_string(md_ij));
                        Vector3 delta_ri = n_ij * 0.5 * (md_ij - d) * res[3];
                        print_debug(delta_ri);
                        Vector3 delta_ri_plus = n_ij * 0.5 * (md_ij - d) * (1 -  res[3]);
                        print_debug(delta_ri_plus);
                        Vector3 delta_rj = n_ij * 0.5 * (d - md_ij) * res[4];
                        print_debug(delta_rj);
                        Vector3 delta_rj_plus = n_ij * 0.5 * (d - md_ij) * (1 -  res[4]);
                        print_debug(delta_rj_plus);

                        m_vertex_velocities.col(i) += delta_ri / delta_time;
                        m_vertex_positions.col(i) += delta_ri;

                        m_vertex_velocities.col(i+1) += delta_ri_plus / delta_time;
                        m_vertex_positions.col(i+1) += delta_ri_plus;

                        other.m_vertex_velocities.col(j) += delta_rj / delta_time;
                        other.m_vertex_positions.col(j) += delta_rj;

                        other.m_vertex_velocities.col(j+1) += delta_rj_plus / delta_time;
                        other.m_vertex_positions.col(j+1) += delta_rj_plus;

                        res = minimum_distance(i, j, m_vertex_positions, other.m_vertex_positions);
                        n_ij = res.head(3);
                        md_ij = n_ij.norm();
                        err_vec[i * (other.m_vertex_positions.cols() - 1) + j] = d - md_ij;
                    } else {
                        err_vec[i * (other.m_vertex_positions.cols() - 1) + j] = 0;
                    }
                } else {
                    err_vec[i * (other.m_vertex_positions.cols() - 1) + j] = 0;
                }
            }
        }
        print_debug("After loop" + std::to_string(k) + "The current error is:" + std::to_string(err_vec.squaredNorm()));

    }

}
