// components/Splitter.cpp
#include "Splitter.h"
#include <iostream> // For warning messages (though none are currently used in this implementation)

/**
 * @brief Constructor for the Splitter component.
 *
 * Initializes the splitter with its name and connected nodes.
 * It represents a conceptual "splitter" or a very low-resistance connection,
 * effectively acting as a near-ideal short circuit. The internal resistance
 * is fixed at R_SPLIT.
 *
 * @param name The unique name of the splitter.
 * @param node1 The name of the first connected node.
 * @param node2 The name of the second connected node.
 * @param tolerance_percent Optional tolerance percentage for component value.
 */
Splitter::Splitter(const std::string& name, const std::string& node1, const std::string& node2,
                   double tolerance_percent)
    : Component(name, node1, node2, tolerance_percent) // Call base class constructor
{
    // R_SPLIT is a const member, so it's initialized directly in the class definition.
    // No specific initialization needed here beyond the base class.
    // A warning could be added here if there were specific conditions for the splitter.
}

/**
 * @brief Applies the stamps of the splitter to the MNA matrix (A) and vector (b).
 *
 * This method treats the splitter as a resistor with a very small, fixed resistance
 * (R_SPLIT) to simulate a near-ideal short circuit. The stamps are applied based
 * on the conductance (G = 1/R_SPLIT).
 *
 * @param A The MNA matrix to which stamps are applied.
 * @param b The MNA right-hand side vector to which stamps are applied.
 * @param x_current_guess The current guess for node voltages and branch currents (not used for this static component).
 * @param prev_solution The solution from the previous time step (not used for this static component).
 * @param time The current simulation time (not used for this static component).
 * @param dt The time step size (not used for this static component).
 */
void Splitter::getStamps(
    Eigen::MatrixXd& A, Eigen::VectorXd& b,
    const Eigen::VectorXd& x_current_guess, const Eigen::VectorXd& prev_solution,
    double time, double dt
) {
    // Get the global indices for the two connected nodes.
    // getNodeIndex is assumed to be a method of the base Component class
    // that correctly maps node names to their corresponding indices in the MNA matrix.
    int idx1 = getNodeIndex(node1);
    int idx2 = getNodeIndex(node2);

    // Calculate the conductance of the splitter.
    // Since R_SPLIT is very small (e.g., 1e-9 Ohm), G_SPLIT will be very large (e.g., 1e9 Siemens).
    double G_SPLIT = 1.0 / R_SPLIT;

    // Apply the stamps to the MNA matrix A.
    // For a conductance G between node1 and node2, the stamps are:
    // A(idx1, idx1) += G
    // A(idx2, idx2) += G
    // A(idx1, idx2) -= G
    // A(idx2, idx1) -= G

    // Ensure node indices are valid before stamping (e.g., not -1 if a node is undefined).
    // This check helps prevent out-of-bounds access if node names are not found.
    if (idx1 != -1) {
        A(idx1, idx1) += G_SPLIT;
    }
    if (idx2 != -1) {
        A(idx2, idx2) += G_SPLIT;
    }
    // Apply cross-stamps only if both nodes are valid and distinct.
    // If idx1 == idx2, it means the component is connected to the same node twice,
    // which is usually an error or a degenerate case.
    if (idx1 != -1 && idx2 != -1 && idx1 != idx2) {
        A(idx1, idx2) -= G_SPLIT;
        A(idx2, idx1) -= G_SPLIT;
    }

    // For a passive, linear component like a resistor (or a very low resistance),
    // there are no contributions to the right-hand side vector 'b' from the component itself,
    // unless there's an independent source associated with it (which is not the case here).
    // Therefore, 'b' remains unchanged by the splitter's stamps.
}

/**
 * @brief Updates the internal state of the splitter.
 *
 * As the splitter is modeled as a static, near-ideal connection (its behavior
 * depends only on the instantaneous voltages across it, not on past states),
 * this method does not perform any state updates. It is provided to satisfy
 * the Component base class interface, ensuring all virtual methods are implemented.
 *
 * @param v_curr The current voltage across the component (not used).
 * @param i_curr The current flowing through the component (not used).
 */
void Splitter::updateState(double v_curr, double i_curr) {
    // No internal state to update for this splitter model.
    // This method is intentionally left empty.
}
