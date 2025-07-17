// components/Resistor.cpp
#include "Resistor.h"
// Assuming utils/Random.h is available if applyTolerance is used.
// If not, you might need to implement applyTolerance directly or remove it.
// #include "../utils/Random.h" // Per applyTolerance
#include <stdexcept>

// Placeholder for applyTolerance if Random.h is not available or desired.
// In a real scenario, this would be a utility function or a member of a base class.
double applyTolerance(double nominal_value, double tolerance_percent) {
    // For simplicity, without Random.h, we return the nominal value.
    // If you need actual tolerance, you'd add random variation here.
    return nominal_value;
}


Resistor::Resistor(const std::string& name, const std::string& node1, const std::string& node2,
                   double resistance_nominal, double tolerance_percent)
    : Component(name, {node1, node2}) {
    resistance = applyTolerance(resistance_nominal, tolerance_percent);
    if (resistance <= 0) {
        throw std::runtime_error("Resistance must be positive for Resistor " + name);
    }
}

/**
 * @brief Applies the stamps of the resistor to the MNA matrix (A) and vector (b).
 *
 * For a resistor, the stamps are constant and linear, based on Ohm's law.
 *
 * @param A The MNA matrix to which stamps are applied.
 * @param b The MNA right-hand side vector to which stamps are applied.
 * @param x_current_guess The current guess for node voltages and branch currents (not used for linear resistors).
 * @param prev_solution The solution from the previous time step (not used for static resistors).
 * @param time The current simulation time (not used for static resistors).
 * @param dt The time step size (not used for static resistors).
 */
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

    // Contributions to the MNA matrix (admittances)
    // If a node is ground (ID 0), its row/column is not stamped.
    if (node1_id != 0) A(node1_id, node1_id) += G;
    if (node2_id != 0) A(node2_id, node2_id) += G;
    if (node1_id != 0 && node2_id != 0) { // Only apply cross-terms if both nodes are not ground
        A(node1_id, node2_id) -= G;
        A(node2_id, node1_id) -= G;
    }
}

/**
 * @brief Updates the internal state of the resistor.
 *
 * For an ideal resistor, there is no internal state that evolves over time.
 * This method is left empty for this model.
 *
 * @param v_curr The current voltage across the component.
 * @param i_curr The current flowing through the component.
 */
void Resistor::updateState(double v_curr, double i_curr) {
    // No state to update for this static linear component.
}
