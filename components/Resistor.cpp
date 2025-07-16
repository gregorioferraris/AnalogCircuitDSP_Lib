// components/Resistor.cpp
#include "Resistor.h"
#include "../utils/Random.h" // Per applyTolerance
#include <stdexcept>

Resistor::Resistor(const std::string& name, const std::string& node1, const std::string& node2,
                   double resistance_nominal, double tolerance_percent)
    : Component(name, {node1, node2}) {
    resistance = applyTolerance(resistance_nominal, tolerance_percent);
    if (resistance <= 0) {
        throw std::runtime_error("Resistance must be positive for Resistor " + name);
    }
}

void Resistor::getStamps(
    Eigen::MatrixXd& A, Eigen::VectorXd& b,
    const Eigen::VectorXd& x_current_guess, const Eigen::VectorXd& prev_solution,
    double time, double dt
) {
    if (nodeIdsVector.size() != 2) {
        throw std::runtime_error("Resistor " + name + " expects 2 node IDs.");
    }
    int node1_id = nodeIdsVector[0];
    int node2_id = nodeIdsVector[1];

    double G = 1.0 / resistance;

    // Contributi alla matrice MNA (ammettenze)
    if (node1_id != 0) A(node1_id, node1_id) += G;
    if (node2_id != 0) A(node2_id, node2_id) += G;
    if (node1_id != 0 && node2_id != 0) {
        A(node1_id, node2_id) -= G;
        A(node2_id, node1_id) -= G;
    }
}
