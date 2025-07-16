// components/JFET.cpp
#include "JFET.h"
#include <iostream> // For error messages
#include <limits>   // For std::numeric_limits

/**
 * @brief Constructor for the JFET component.
 *
 * Initializes an N-channel JFET with its name, three connected nodes (Gate, Drain, Source),
 * and key parameters.
 *
 * @param name The unique name of the JFET.
 * @param node_gate The name of the Gate node.
 * @param node_drain The name of the Drain node.
 * @param node_source The name of the Source node.
 * @param Idss The drain current with VGS=0 (Amps). Must be > 0.
 * @param Vp The pinch-off voltage (Volts). Must be < 0 for N-channel.
 * @param tolerance_percent Optional tolerance percentage.
 */
JFET::JFET(const std::string& name,
           const std::string& node_gate, const std::string& node_drain, const std::string& node_source,
           double Idss, double Vp, double tolerance_percent)
    // The base Component constructor takes node1 and node2.
    // For a 3-terminal device, we can map node1 to Drain and node2 to Source,
    // and handle node_gate separately. The base Component's node1 and node2
    // are not strictly used for the JFET's internal connections, but for its
    // overall presence in the circuit.
    : Component(name, node_drain, node_source, tolerance_percent), // node1=drain, node2=source
      node_gate(node_gate), node_drain(node_drain), node_source(node_source),
      Idss(Idss), Vp(Vp),
      prev_VGS(0.0) // Initialize state variable
{
    // A JFET model typically does not require additional auxiliary variables
    // beyond the node voltages when using companion models for non-linearities.
    setNumAuxiliaryVariables(0);

    // Validate input parameters
    if (Idss <= 0.0) {
        std::cerr << "Warning: JFET " << name << " has a non-positive Idss. Setting to 1e-3A." << std::endl;
        this->Idss = 1e-3; // Default typical value
    }
    // For N-channel JFET, Vp should be negative.
    if (Vp >= 0.0) {
        std::cerr << "Warning: JFET " << name << " has a non-negative Vp (pinch-off voltage). Setting to -2.0V (assuming N-channel)." << std::endl;
        this->Vp = -2.0; // Default typical value for N-channel
    }

    std::cout << "JFET " << name << " initialized with Idss=" << Idss << "A, Vp=" << Vp << "V." << std::endl;
}

/**
 * @brief Helper function to calculate the drain current (ID) based on VGS.
 *
 * Implements the square-law characteristic: ID = Idss * (1 - VGS / Vp)^2.
 * Handles the cutoff region (VGS >= Vp) where ID = 0.
 *
 * @param VGS The Gate-Source voltage.
 * @return The calculated drain current.
 */
double JFET::calculateDrainCurrent(double VGS) {
    if (VGS >= Vp) { // Cutoff region
        return 0.0;
    }
    // Square-law model for saturation region
    double term = (1.0 - VGS / Vp);
    return Idss * term * term;
}

/**
 * @brief Helper function to calculate the transconductance (gm) based on VGS.
 *
 * gm = d(ID)/d(VGS) = -2 * Idss / Vp * (1 - VGS / Vp).
 * Handles the cutoff region (VGS >= Vp) where gm = 0.
 *
 * @param VGS The Gate-Source voltage.
 * @return The calculated transconductance.
 */
double JFET::calculateTransconductance(double VGS) {
    if (VGS >= Vp) { // Cutoff region
        return 0.0;
    }
    // Derivative of the square-law model
    return -2.0 * Idss / Vp * (1.0 - VGS / Vp);
}

/**
 * @brief Applies the stamps of the JFET model to the MNA matrix (A) and vector (b).
 *
 * This method implements the linearized companion model for the JFET's
 * voltage-controlled drain current. It calculates the transconductance (gm)
 * and an equivalent current source based on the previous time step's VGS.
 * The current flows from Drain to Source.
 *
 * @param A The MNA matrix to which stamps are applied.
 * @param b The MNA right-hand side vector to which stamps are applied.
 * @param x_current_guess The current guess for node voltages and branch currents (used for updating state).
 * @param prev_solution The solution from the previous time step (used for companion models).
 * @param time The current simulation time (not directly used for stamping).
 * @param dt The time step size (not directly used for this static model, but implicitly affects prev_solution).
 */
void JFET::getStamps(
    Eigen::MatrixXd& A, Eigen::VectorXd& b,
    const Eigen::VectorXd& x_current_guess, const Eigen::VectorXd& prev_solution,
    double time, double dt
) {
    // Get the global indices for the JFET nodes.
    // Assumes getNodeIndex is available and correctly maps node names to matrix indices.
    int idx_G = getNodeIndex(node_gate);
    int idx_D = getNodeIndex(node_drain);
    int idx_S = getNodeIndex(node_source);

    // Check for valid indices
    if (idx_G == -1 || idx_D == -1 || idx_S == -1) {
        std::cerr << "Error: Invalid node index for JFET " << name << std::endl;
        return;
    }

    // Calculate gm and the equivalent current source based on prev_VGS
    double gm = calculateTransconductance(prev_VGS);
    double Id_prev = calculateDrainCurrent(prev_VGS);

    // The linearized current is: I_D = gm * VGS + I_eq
    // where I_eq = Id_prev - gm * prev_VGS
    // And VGS = V(idx_G) - V(idx_S)

    // Stamp the gm * VGS term: gm * (V(idx_G) - V(idx_S))
    // This current flows from Drain (idx_D) to Source (idx_S).
    // It's subtracted from the Drain node equation and added to the Source node equation.
    A(idx_D, idx_G) -= gm; // Contribution from V(idx_G) to I_D
    A(idx_D, idx_S) += gm; // Contribution from V(idx_S) to I_D

    A(idx_S, idx_G) += gm; // Contribution from V(idx_G) to I_S (which is -I_D)
    A(idx_S, idx_S) -= gm; // Contribution from V(idx_S) to I_S

    // Stamp the constant current source term: I_eq = Id_prev - gm * prev_VGS
    double I_eq = Id_prev - gm * prev_VGS;

    // This current flows from Drain (idx_D) to Source (idx_S).
    // It's added to the Drain node equation and subtracted from the Source node equation.
    b(idx_D) -= I_eq; // Current flows OUT of Drain
    b(idx_S) += I_eq; // Current flows INTO Source
}

/**
 * @brief Updates the internal state variable (prev_VGS) for the next time step.
 *
 * This method is called after the MNA system is solved for the current time step.
 * It extracts the current Gate-Source voltage from the solution vector and stores it
 * for use in the companion model of the next time step.
 *
 * @param v_curr The current voltage across the component (not used directly for JFET state).
 * @param i_curr The current flowing through the component (not used directly).
 */
void JFET::updateState(double v_curr, double i_curr) {
    // The `v_curr` and `i_curr` parameters here are for the component as a whole,
    // which is not directly meaningful for a 3-terminal device like a JFET.
    // Instead, we need to access the solved node voltages from the full solution vector.
    // This typically happens in the main Circuit loop after A*x=b is solved.
    // It is assumed that `x_current_guess` (which is the full solution vector `x` from the current time step)
    // is available and contains the latest node voltages.

    // Get the global indices for the JFET nodes.
    int idx_G = getNodeIndex(node_gate);
    int idx_S = getNodeIndex(node_source);

    // Check if `x_current_guess` is large enough to contain the node voltages.
    // This check is a safeguard; proper integration with the Circuit solver is crucial.
    if (idx_G < 0 || idx_G >= x_current_guess.size() ||
        idx_S < 0 || idx_S >= x_current_guess.size()) {
        std::cerr << "Error: Solution vector size mismatch or invalid node index for JFET " << name << " in updateState." << std::endl;
        return;
    }

    double V_G_curr = x_current_guess(idx_G);
    double V_S_curr = x_current_guess(idx_S);

    // Update the previous VGS for the next time step
    prev_VGS = V_G_curr - V_S_curr;
}
