#include <DiscreteElasticRod.hpp>
#include <common.h>

static VectorX minimum_distance(const int i, const int j, const Matrix3X &vertices) {
    Vector3 xi = Vertices.col(i);
    Vector3 xi_plus = Vertices.col(i+1);
    Vector3 xj = Vertices.col(j);
    Vector3 xj_plus = Vertices.col(j+1);

    Vector3 ei = xi_plus - xi;
    Vector3 ej = xj_plus - xj;

    Vector3 n = ei.cross(ej);
    if(n.squaredNorm < 1e-8) {
        return Vector3::Zero();
    }

    Float t1 = clamp((ej.cross(n).dot(xj - xi)) / n.squaredNorm(), 0, 1);
    Float t2 = clamp((ei.cross(n).dot(xj - xi)) / n.squaredNorm(), 0, 1);

    Vector3 p1 = xi + t1 * ei;
    Vector3 p2 = xj + t2 * ej;

    VectorX res(5);
    res.head(3) = (p2 - p1);
    res[3] = t1;
    res[4] = t2;

    return res;
}

void DiscreteElasticRod::calculateContactForces(const Float d) {



    for (int i = 0; i < m_n; i++) {
        for (int j = 0; j < m_n; j++) {
            if(abs(i - j) > 1) {
                VectorX res = minimum_distance(i, j, m_vertices);
                Vector3 n_ij = res.head(3);
                Float md_ij = n_ij.norm();
                if (md_ij < d) {
                    for (int k = 0; k < max_iter && md_ij < d; k++) {
                        Vector delta_ri = n_ij * 0.5 * (md_ij - d) * res[3];
                        Vector delta_ri_plus = n_ij * 0.5 * (md_ij - d) * (1 -  res[3]);
                        Vector delta_rj = n_ij * 0.5 * (d - md_ij) * res[4];
                        Vector delta_rj_plus = n_ij * 0.5 * (d - md_ij) * (1 -  res[4]);

                        m_vertices.col(i) += delta_ri;
                        m_vertices.col(i+1) += delta_ri_plus;
                        m_vertices.col(j) += delta_rj;
                        m_vertices.col(j+1) += delta_rj_plus;
                        res = minimum_distance(i, j, m_vertices);
                        n_ij = res.head(3);
                        md_ij = n_ij.norm();
                    }
                }
            }
        }
    }

}
