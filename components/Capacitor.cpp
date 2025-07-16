// components/Capacitor.cpp
#include "Capacitor.h"
#include "../utils/Random.h" // Per applyTolerance
#include <stdexcept>

Capacitor::Capacitor(const std::string& name, const std::string& node1, const std::string& node2,
                     double capacitance_nominal, double tolerance_percent)
    : Component(name, {node1, node2}), v_prev(0.0), i_prev(0.0) {
    capacitance = applyTolerance(capacitance_nominal, tolerance_percent);
    if (capacitance <= 0) {
        throw std::runtime_error("Capacitance must be positive for Capacitor " + name);
    }
}

void Capacitor::getStamps(
    Eigen::MatrixXd& A, Eigen::VectorXd& b,
    const Eigen::VectorXd& x_current_guess, const Eigen::VectorXd& prev_solution,
    double time, double dt
) {
    if (nodeIdsVector.size() != 2) {
        throw std::runtime_error("Capacitor " + name + " expects 2 node IDs.");
    }
    int node1_id = nodeIdsVector[0];
    int node2_id = nodeIdsVector[1];

    double G_eq = 2.0 * capacitance / dt;

    // Contributi alla matrice MNA (parte dipendente dalla tensione attuale)
    if (node1_id != 0) A(node1_id, node1_id) += G_eq;
    if (node2_id != 0) A(node2_id, node2_id) += G_eq;
    if (node1_id != 0 && node2_id != 0) {
        A(node1_id, node2_id) -= G_eq;
        A(node2_id, node1_id) -= G_eq;
    }

    // Contributi al vettore RHS (parte dipendente dallo stato precedente)
    double V_C_prev = prev_solution[node1_id] - prev_solution[node2_id];
    double i_eq = G_eq * V_C_prev + i_prev; // i_prev è la corrente al passo precedente

    if (node1_id != 0) b(node1_id) -= i_eq;
    if (node2_id != 0) b(node2_id) += i_eq;
}

void Capacitor::updateState(double v_curr, double i_curr) {
    v_prev = v_curr;
    i_prev = i_curr;
}
