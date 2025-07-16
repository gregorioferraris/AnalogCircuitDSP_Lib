// components/Inductor.cpp
#include "Inductor.h"
#include "../utils/Random.h" // Per applyTolerance
#include <stdexcept>

Inductor::Inductor(const std::string& name, const std::string& node1, const std::string& node2,
                   double inductance_nominal, double tolerance_percent)
    : Component(name, {node1, node2}), i_prev(0.0), v_prev(0.0) {
    L = applyTolerance(inductance_nominal, tolerance_percent);
    if (L <= 0) {
        throw std::runtime_error("Inductance must be positive for Inductor " + name);
    }
}

void Inductor::getStamps(
    Eigen::MatrixXd& A, Eigen::VectorXd& b,
    const Eigen::VectorXd& x_current_guess, const Eigen::VectorXd& prev_solution,
    double time, double dt
) {
    if (nodeIdsVector.size() != 2) {
        throw std::runtime_error("Inductor " + name + " expects 2 node IDs.");
    }
    int node1_id = nodeIdsVector[0];
    int node2_id = nodeIdsVector[1];

    double G_eq = dt / (2.0 * L);

    // Contributi alla matrice MNA
    if (node1_id != 0) A(node1_id, node1_id) += G_eq;
    if (node2_id != 0) A(node2_id, node2_id) += G_eq;
    if (node1_id != 0 && node2_id != 0) {
        A(node1_id, node2_id) -= G_eq;
        A(node2_id, node1_id) -= G_eq;
    }

    // Contributi al vettore RHS (parte dipendente dallo stato precedente)
    double V_eq = i_prev * (2.0 * L / dt) + v_prev;
    double i_eq = G_eq * V_eq;

    if (node1_id != 0) b(node1_id) -= i_eq;
    if (node2_id != 0) b(node2_id) += i_eq;
}

void Inductor::updateState(double v_curr, double i_curr) {
    v_prev = v_curr;
    i_prev = i_curr;
}
