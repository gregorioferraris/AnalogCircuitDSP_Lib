// components/CloseBoxCabinet.cpp
#include "CloseBoxCabinet.h"
#include <iostream> // For output messages

/**
 * @brief Constructor for the CloseBoxCabinet component.
 *
 * Initializes a Close Box Cabinet with its name and two connected nodes.
 *
 * @param name The unique name of the cabinet.
 * @param node1 The name of the first conceptual node associated with the cabinet.
 * @param node2 The name of the second conceptual node associated with the cabinet.
 * @param tolerance_percent Optional tolerance percentage (inherited from Component).
 */
CloseBoxCabinet::CloseBoxCabinet(const std::string& name,
                                 const std::string& node1, const std::string& node2,
                                 double tolerance_percent)
    : Component(name, node1, node2, tolerance_percent)
{
    // A CloseBoxCabinet typically does not add auxiliary variables unless
    // it models complex internal electrical behavior (e.g., internal shielding currents).
    setNumAuxiliaryVariables(0);

    std::cout << "CloseBoxCabinet " << name << " initialized. This component is primarily for organizational purposes." << std::endl;
}

/**
 * @brief Applies the stamps of the CloseBoxCabinet to the MNA matrix (A) and vector (b).
 *
 * For a basic CloseBoxCabinet, this method does not add any stamps to the MNA matrix
 * or vector. It serves as a placeholder or for conceptual grouping of components.
 *
 * If specific electrical effects (e.g., parasitic capacitance to ground,
 * shielding effects, or defining a common internal ground plane) need to be modeled,
 * this method would be extended to include the appropriate stamps.
 *
 * @param A The MNA matrix to which stamps are applied.
 * @param b The MNA right-hand side vector to which stamps are applied.
 * @param x_current_guess The current guess for node voltages and branch currents (not used).
 * @param prev_solution The solution from the previous time step (not used).
 * @param time The current simulation time (not used).
 * @param dt The time step size (not used).
 */
void CloseBoxCabinet::getStamps(
    Eigen::MatrixXd& A, Eigen::VectorXd& b,
    const Eigen::VectorXd& x_current_guess, const Eigen::VectorXd& prev_solution,
    double time, double dt
) {
    // This component typically does not add stamps to the MNA matrix.
    // It's mainly for organizational purposes or if complex parasitic/shielding
    // effects are modeled, which would require specific implementations here.

    // Example of what *could* be done if it represented a direct connection
    // to ground for one of its nodes (e.g., node2 is internally grounded):
    // int idx2 = getNodeIndex(node2);
    // if (idx2 != -1) {
    //     A(idx2, idx2) += 1.0; // Connects node2 to ground (assuming ground is node 0 or implicit)
    // }
    // However, for a generic "cabinet", this is usually not the default behavior.

    // No operations on A or b for a basic CloseBoxCabinet.
}
