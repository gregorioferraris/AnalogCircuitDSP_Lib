// components/OpAmp.cpp
#include "OpAmp.h"
#include <iostream> // For debug output and warnings

/**
 * @brief Constructor for the OpAmp component.
 *
 * Initializes an ideal operational amplifier with its name and three connection nodes.
 *
 * @param name The unique name of the OpAmp.
 * @param positive_input_node The name of the non-inverting (+) input node.
 * @param negative_input_node The name of the inverting (-) input node.
 * @param output_node The name of the output node.
 * @param tolerance_percent Optional tolerance percentage (inherited from Component).
 */
OpAmp::OpAmp(const std::string& name,
             const std::string& positive_input_node,
             const std::string& negative_input_node,
             const std::string& output_node,
             double tolerance_percent)
    : Component(name, positive_input_node, negative_input_node, tolerance_percent), // Pass two nodes to base Component (arbitrary for OpAmp, but needed for constructor)
      pos_input_node_name(positive_input_node),
      neg_input_node_name(negative_input_node),
      out_node_name(output_node)
{
    // An ideal OpAmp introduces one auxiliary variable for its output current.
    setNumAuxiliaryVariables(1);

    // Add the output node to the list of nodes managed by the component.
    // This is crucial for getNodeIndex to work correctly for the output node.
    addNode(output_node);

    std::cout << "OpAmp " << name << " initialized. Inputs: " << positive_input_node
              << " (+), " << negative_input_node << " (-). Output: " << output_node << "." << std::endl;
}

/**
 * @brief Applies the "stamps" of the OpAmp to the MNA matrix (A) and vector (b).
 *
 * This method implements the ideal OpAmp model.
 *
 * The ideal OpAmp model enforces two conditions:
 * 1. No current flows into the input terminals (infinite input impedance).
 * This means the KCL equations for the positive and negative input nodes are not directly affected
 * by the OpAmp itself, other than through external connections.
 * 2. The voltage difference between the input terminals is zero (infinite gain).
 * This is enforced by adding a constraint equation to the MNA system.
 * V_pos - V_neg = 0
 *
 * An auxiliary variable is introduced for the output current (I_out) flowing out of the output terminal.
 *
 * MNA stamps:
 * 1. KCL at output node (idx_out_1based):
 * The output current I_out flows OUT of the output node.
 * A(idx_out_1based, idx_opamp_current_1based) += 1.0;
 *
 * 2. Constraint equation (auxiliary row idx_opamp_current_1based):
 * V_pos - V_neg = 0
 * A(idx_opamp_current_1based, idx_p_1based) += 1.0;
 * A(idx_opamp_current_1based, idx_n_1based) += -1.0;
 * b(idx_opamp_current_1based) += 0.0;
 *
 * @param A The MNA matrix to which the stamps are applied.
 * @param b The MNA right-hand side vector to which the stamps are applied.
 * @param x_current_guess The current guess for node voltages and branch currents (not directly used for this linear model).
 * @param prev_solution The solution from the previous time step (not used for this static component).
 * @param time The current simulation time (not used for this static component).
 * @param dt The time step size (not used for this static component).
 */
void OpAmp::getStamps(
    Eigen::MatrixXd& A, Eigen::VectorXd& b,
    const Eigen::VectorXd& x_current_guess, const Eigen::VectorXd& prev_solution,
    double time, double dt
) {
    // Get the global (1-based) indices for the input and output nodes.
    // Assuming getNodeIndex(node_name) returns the 1-based index (0 for ground).
    int idx_p_1based = getNodeIndex(pos_input_node_name);
    int idx_n_1based = getNodeIndex(neg_input_node_name);
    int idx_out_1based = getNodeIndex(out_node_name);

    // Get the global (1-based) index for the auxiliary current variable introduced by this component.
    // Assuming getAuxiliaryVariableStartIndex() returns the 1-based index of the first (and only)
    // auxiliary variable for this OpAmp component.
    int idx_opamp_current_1based = getAuxiliaryVariableStartIndex();

    // --- Apply MNA stamps ---

    // 1. KCL equation for the output node (idx_out_1based):
    // The output current I_out flows OUT of the output node. So, its coefficient is +1.
    // A(row_idx, col_idx) += value
    A(idx_out_1based, idx_opamp_current_1based) += 1.0;

    // 2. Constraint equation for the OpAmp (auxiliary row idx_opamp_current_1based):
    // Enforce V_pos - V_neg = 0
    // Coefficient for V_pos is +1
    A(idx_opamp_current_1based, idx_p_1based) += 1.0;
    // Coefficient for V_neg is -1
    A(idx_opamp_current_1based, idx_n_1based) += -1.0;

    // The right-hand side for this constraint equation is 0.
    b(idx_opamp_current_1based) += 0.0; // Explicitly adding 0.0 for clarity, though it's often implicitly 0.
}

/**
 * @brief Updates the internal state of the OpAmp.
 *
 * For an ideal OpAmp, there is no internal state to update
 * based on previous time steps. This method is empty.
 *
 * @param v_curr The current voltage across the component (not used).
 * @param i_curr The current flowing through the component (not used).
 */
void OpAmp::updateState(double v_curr, double i_curr) {
    // No state to update for this static non-linear model.
}
