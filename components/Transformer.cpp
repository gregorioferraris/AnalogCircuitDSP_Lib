// components/Transformer.cpp
#include "Transformer.h"
#include <iostream> // For error messages

/**
 * @brief Constructor for the Transformer component.
 *
 * Initializes an ideal transformer with its name, four connected nodes,
 * and its turns ratio.
 *
 * @param name The unique name of the transformer.
 * @param node1_p The name of the first node on the primary side (positive terminal).
 * @param node2_p The name of the second node on the primary side (negative terminal).
 * @param node1_s The name of the first node on the secondary side (positive terminal).
 * @param node2_s The name of the second node on the secondary side (negative terminal).
 * @param turns_ratio The turns ratio (Np/Ns) of the transformer. Must be non-zero.
 * @param tolerance_percent Optional tolerance percentage for component value.
 */
Transformer::Transformer(const std::string& name, const std::string& node1_p, const std::string& node2_p,
                         const std::string& node1_s, const std::string& node2_s,
                         double turns_ratio, double tolerance_percent)
    // Pass primary nodes to the base Component constructor.
    // The Component base class is assumed to handle node1 and node2.
    // The secondary nodes are managed directly by the Transformer class.
    : Component(name, node1_p, node2_p, tolerance_percent),
      node1_s(node1_s), node2_s(node2_s), turns_ratio(turns_ratio)
{
    // An ideal transformer requires two auxiliary variables for its primary and secondary currents.
    // We assume the base Component class has a mechanism (e.g., setNumAuxiliaryVariables)
    // to request multiple auxiliary variables, and the simulator will assign contiguous indices.
    setNumAuxiliaryVariables(2);

    if (turns_ratio == 0.0) {
        std::cerr << "Warning: Transformer " << name << " has a turns ratio of zero. This can lead to division by zero." << std::endl;
    }
}

/**
 * @brief Applies the stamps of the ideal transformer to the MNA matrix (A) and vector (b).
 *
 * This method introduces two auxiliary variables for the primary (Ip) and secondary (Is)
 * currents. It stamps the voltage and current relationships of an ideal transformer.
 *
 * @param A The MNA matrix to which stamps are applied.
 * @param b The MNA right-hand side vector to which stamps are applied.
 * @param x_current_guess The current guess for node voltages and branch currents (not used for this static component).
 * @param prev_solution The solution from the previous time step (not used for this static component).
 * @param time The current simulation time (not used for this static component).
 * @param dt The time step size (not used for this static component).
 */
void Transformer::getStamps(
    Eigen::MatrixXd& A, Eigen::VectorXd& b,
    const Eigen::VectorXd& x_current_guess, const Eigen::VectorXd& prev_solution,
    double time, double dt
) {
    // Get the global indices for all four connected nodes.
    // getNodeIndex is assumed to be a method of the base Component class (or accessible via Circuit)
    // that correctly maps node names to their corresponding indices in the MNA matrix.
    int idx_p1 = getNodeIndex(node1); // Primary positive node (inherited from Component::node1)
    int idx_p2 = getNodeIndex(node2); // Primary negative node (inherited from Component::node2)
    int idx_s1 = getNodeIndex(node1_s); // Secondary positive node
    int idx_s2 = getNodeIndex(node2_s); // Secondary negative node

    // Get the starting index for the auxiliary variables.
    // The first auxiliary variable (for Ip) will be at aux_start_idx.
    // The second auxiliary variable (for Is) will be at aux_start_idx + 1.
    int aux_start_idx = getAuxiliaryVariableStartIndex();
    int aux_idx_Ip = aux_start_idx;
    int aux_idx_Is = aux_start_idx + 1;

    // Check for valid indices. If any index is -1, it means a node or auxiliary variable
    // was not properly registered/found in the global mapping.
    if (idx_p1 == -1 || idx_p2 == -1 || idx_s1 == -1 || idx_s2 == -1 || aux_start_idx == -1) {
        std::cerr << "Error: Invalid node or auxiliary variable index for Transformer " << name << std::endl;
        return;
    }

    // Ensure turns_ratio is not zero to avoid division by zero.
    if (turns_ratio == 0.0) {
        // This case should ideally be caught in the constructor or handled gracefully.
        // For stamping, we can effectively treat it as an open circuit on the secondary
        // or prevent stamping if it's truly problematic. For now, just return.
        std::cerr << "Error: Transformer " << name << " has a turns ratio of zero, cannot stamp." << std::endl;
        return;
    }

    // --- Apply stamps for the ideal transformer equations ---

    // Equation 1: V_primary - a * V_secondary = 0
    // (V(node1_p) - V(node2_p)) - turns_ratio * (V(node1_s) - V(node2_s)) = 0
    // This equation is placed in the row corresponding to the primary current auxiliary variable (Ip).
    A(aux_idx_Ip, idx_p1) += 1.0;
    A(aux_idx_Ip, idx_p2) -= 1.0;
    A(aux_idx_Ip, idx_s1) -= turns_ratio;
    A(aux_idx_Ip, idx_s2) += turns_ratio;
    // b(aux_idx_Ip) remains 0.0 as it's a homogeneous equation.

    // Equation 2: Ip + (1/a) * Is = 0
    // This equation is placed in the row corresponding to the secondary current auxiliary variable (Is).
    A(aux_idx_Is, aux_idx_Ip) += 1.0;
    A(aux_idx_Is, aux_idx_Is) += 1.0 / turns_ratio;
    // b(aux_idx_Is) remains 0.0 as it's a homogeneous equation.

    // --- Apply current contributions to node equations ---

    // Primary side current: Ip enters node1_p and leaves node2_p
    A(idx_p1, aux_idx_Ip) += 1.0;
    A(idx_p2, aux_idx_Ip) -= 1.0;

    // Secondary side current: Is enters node1_s and leaves node2_s
    A(idx_s1, aux_idx_Is) += 1.0;
    A(idx_s2, aux_idx_Is) -= 1.0;
}

/**
 * @brief Updates the internal state of the transformer.
 *
 * As an ideal transformer is a static, lossless component (its behavior
 * depends only on instantaneous voltages and currents, not on past states),
 * this method does not perform any state updates. It is provided to satisfy
 * the Component base class interface, ensuring all virtual methods are implemented.
 *
 * @param v_curr The current voltage across the component (not used).
 * @param i_curr The current flowing through the component (not used).
 */
void Transformer::updateState(double v_curr, double i_curr) {
    // No internal state to update for an ideal transformer model.
    // This method is intentionally left empty.
}
