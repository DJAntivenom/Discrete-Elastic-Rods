#include <DiscreteElasticRod.hpp>
#include <Optimization.h>

#include <stdexcept>

bool DiscreteElasticRod::getConstraints(const Optimization::VectorXf &q_r_x, double &energy) const {

    Eigen::VectorXf m_mass = Eigen::MatrixXf::Ones(m_vertex_positions.size()) / 3.f;

    Eigen::MatrixXf M = Eigen::MatrixXf::Zero(3 * m_n + 12, 3 * m_n + 12);
    M.block<3, 3>(0, 0) = 4.f * Eigen::MatrixXf::Identity(3, 3);
    M.block<3, 3>(3, 3) = m_mass.sum() * Eigen::MatrixXf::Identity(3, 3);

    M.diagonal().tail(3 * m_n + 6) = m_mass; 

    Eigen::Quaternionf q = Eigen::Quaternionf(1, 0, 0, 0);

    //TODO: get q_deriv apparently there is also a q.derived() function? but that can't do what i think it does, can it?
    Eigen::Quaternionf q_deriv = Eigen::Quaternionf(1, 0, 0, 0); 

    Eigen::Vector3f r = Eigen::Vector3f::Zero();
    //TODO: and how are we supposed to get the derivative of r in time either?
    Eigen::Vector3f r_deriv = Eigen::Vector3f(1, 0, 0);

    Eigen::VectorXf y = Eigen::VectorXf::Zero(3 * m_n + 12);
    y.segment<3>(0) = (q.inverse() * q_deriv).coeffsScalarFirst().segment<3>(1);
    y.segment<3>(3) = r_deriv;
    y.tail(3 * m_n + 6) = m_vertex_velocities;

    Eigen::VectorXf C = Eigen::VectorXf::Zero(m_n + 7);
    for (int i = 0; i <= m_n; i++) {
        C(i) = (m_vertex_positions.block(3, i+1, 3, 1) - m_vertex_positions.block(3, i, 3, 1)).squaredNorm() - m_edge_length(i) * m_edge_length(i);
    }
    C(m_n + 1) = q.squaredNorm() - 1.f;
    //TODO: do we need to save rest positions for all vertices?
    C.segment<3>(m_n + 2) = q._transformVector(m_vertex_positions.segment<3>(0)) + r - m_vertex_positions.segment<3>(0); 
    C.segment<3>(m_n + 5) = q._transformVector(m_vertex_positions.segment<3>(1)) + r - m_vertex_positions.segment<3>(1);
    
    //TODO: how do we get the values for lambda? is that passed in with the optimization vector?
    Eigen::VectorXf lambda = Eigen::VectorXf::Zero(m_n + 7);

    //this is assuming that the output value is a scalar (using y^TMy instead of the yMy^T written in the paper), because otherwise the dimensions of this calculation don't work
    energy = y.transpose() * M * y - C.dot(lambda);


    return true;
}

bool DiscreteElasticRod::getConstraintGradient(const Optimization::VectorXf &q_r_x, Optimization::VectorXf &gradient) const {

    //TODO: Figure out how to take the derivative of y^TMy - Cl
    //contribution from inextensibility constraints
    Eigen::VectorXf X_constraints = Eigen::VectorXf::Zero(3 * (m_n + 2));
    for (int i = 0; i <= m_n; i++) {
        X_constraints.segment<3>(3 * i) = -(m_vertex_positions.block(3, i+1, 3, 1) - m_vertex_positions.block(3, i, 3, 1));
    }
    for (int i = 1; i <= m_n + 1; i++) {
        X_constraints.segment<3>(3 * i) += (m_vertex_positions.block(3, i, 3, 1) - m_vertex_positions.block(3, i-1, 3, 1));
    }
    //contribution from rigid-body coupling constraints
    X_constraints.segment<6>(0) -= Eigen::MatrixXf::Ones(6);

    //contribution from generalized kinetic energy
    //X_constraints += m_vertex_mass;

    Eigen::Quaternionf q = Eigen::Quaternionf(1, 0, 0, 0);

    Eigen::VectorXf C = Eigen::VectorXf::Zero(m_n + 7);
    for (int i = 0; i <= m_n; i++) {
        C(i) = (m_vertex_positions.block(3, i+1, 3, 1) - m_vertex_positions.block(3, i, 3, 1)).squaredNorm() - m_edge_length(i) * m_edge_length(i);
    }
    C(m_n + 1) = q.squaredNorm() - 1.f;
    //TODO: do we need to save rest positions for all vertices?
    C.segment<3>(m_n + 2) = q._transformVector(m_vertex_positions.segment<3>(0)) + r - m_vertex_positions.segment<3>(0); 
    C.segment<3>(m_n + 5) = q._transformVector(m_vertex_positions.segment<3>(1)) + r - m_vertex_positions.segment<3>(1);
}
    

bool DiscreteElasticRod::getConstraintHessian(const Optimization::VectorXf &q_r_x, Optimization::TripletListF &hessian) const {

    return false;
}