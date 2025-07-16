// components/OpAmp.h
#ifndef OPAMP_H
#define OPAMP_H

#include "Component.h" // Include the base Component class
#include <string>

/**
 * @brief Represents an ideal Operational Amplifier (OpAmp) component in a circuit simulation.
 *
 * This component models an ideal OpAmp with the following properties:
 * 1. Infinite input impedance: No current flows into the non-inverting (+) and inverting (-) input terminals.
 * 2. Infinite open-loop gain: The voltage difference between the input terminals is zero (V+ = V-).
 *
 * It has three electrical terminals: a non-inverting input node (+), an inverting input node (-),
 * and an output node (out).
 *
 * This component introduces one auxiliary variable to represent the output current flowing
 * out of the output terminal.
 */
class OpAmp : public Component {
public:
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
    OpAmp(const std::string& name,
          const std::string& positive_input_node,
          const std::string& negative_input_node,
          const std::string& output_node,
          double tolerance_percent = 0.0);

    /**
     * @brief Applies the "stamps" of the OpAmp to the MNA matrix (A) and vector (b).
     *
     * This method implements the ideal OpAmp model by enforcing the voltage constraint
     * (V+ = V-) and introducing an auxiliary variable for the output current.
     *
     * @param A The MNA matrix to which the stamps are applied.
     * @param b The MNA right-hand side vector to which the stamps are applied.
     * @param x_current_guess The current guess for node voltages and branch currents (not directly used for this linear model).
     * @param prev_solution The solution from the previous time step (not used for this static component).
     * @param time The current simulation time (not used for this static component).
     * @param dt The time step size (not used for this static component).
     */
    void getStamps(
        Eigen::MatrixXd& A, Eigen::VectorXd& b,
        const Eigen::VectorXd& x_current_guess, const Eigen::VectorXd& prev_solution,
        double time, double dt
    ) override;

private:
    // Node names for clarity (already stored in the base Component class, but useful as local copies)
    std::string pos_input_node_name;
    std::string neg_input_node_name;
    std::string out_node_name;
};

#endif // OPAMP_H
