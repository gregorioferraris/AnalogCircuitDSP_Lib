// components/RibbonMicrophone.cpp
#include "RibbonMicrophone.h"
#include <iostream> // For debug output and warnings
#include <cmath> // For std::abs

/**
 * @brief Constructor for the RibbonMicrophone component.
 *
 * Initializes a ribbon microphone with its name, two output nodes,
 * output resistance, output inductance, and a function that provides
 * the instantaneous sound pressure as a function of time.
 *
 * @param name The unique name of the microphone.
 * @param positive_node The name of the positive output node.
 * @param negative_node The name of the negative output node.
 * @param output_resistance The series output resistance of the microphone (Ohms).
 * @param output_inductance The series output inductance of the microphone (Henry).
 * @param sensitivity The sensitivity of the microphone (Volts/Pascal).
 * @param sound_pressure_function A std::function that takes the current time (double)
 * and returns the instantaneous sound pressure (Pascal).
 * @param tolerance_percent Optional tolerance percentage (inherited from Component).
 */
RibbonMicrophone::RibbonMicrophone(const std::string& name,
                                   const std::string& positive_node,
                                   const std::string& negative_node,
                                   double output_resistance,
                                   double output_inductance,
                                   double sensitivity,
                                   std::function<double(double time)> sound_pressure_function,
                                   double tolerance_percent)
    : Component(name, positive_node, negative_node, tolerance_percent), // Pass the two main nodes to the base Component class
      R_out(output_resistance), L_out(output_inductance),
      sensitivity(sensitivity),
      V_sound_pressure_func(sound_pressure_function),
      pos_node_name(positive_node), neg_node_name(negative_node)
{
    // A ribbon microphone, modeled as a voltage source with series impedance,
    // introduces one auxiliary variable for its branch current in MNA.
    setNumAuxiliaryVariables(1);

    // Basic validation for parameters
    if (R_out < 0 || L_out < 0 || sensitivity < 0) {
        std::cerr << "Warning: RibbonMicrophone " << name << " has negative output resistance (" << R_out << " Ohm), inductance (" << L_out << " Henry) or sensitivity (" << sensitivity << " V/Pa). This might lead to unstable simulations." << std::endl;
    }
    std::cout << "RibbonMicrophone " << name << " initialized. Output between " << positive_node << " and " << negative_node << "." << std::endl;
}

/**
 * @brief Applies the "stamps" of the RibbonMicrophone to the MNA matrix (A) and vector (b).
 *
 * This method implements the companion model for the series combination of a
 * time-varying voltage source, a resistor, and an inductor.
 *
 * The inductor L is modeled using the Backward Euler companion model:
 * V_L(t) = (L/dt) * I_L(t) - (L/dt) * I_L(t-dt)
 *
 * The total branch equation (from positive node to negative node, with current I_mic):
 * V_pos - V_neg = V_mic_generated(t) + I_mic(t) * R_out + V_L(t)
 * Substitute V_L(t):
 * V_pos - V_neg = V_mic_generated(t) + I_mic(t) * R_out + (L_out/dt) * I_mic(t) - (L_out/dt) * I_mic(t-dt)
 *
 * Rearrange into MNA form (terms with unknowns on the left, known terms on the right):
 * V_pos - V_neg - (R_out + L_out/dt) * I_mic(t) = V_mic_generated(t) - (L_out/dt) * I_mic(t-dt)
 *
 * Let R_total_eff = R_out + L_out/dt
 * Let V_eff_source = V_mic_generated(t) - (L_out/dt) * I_mic(t-dt)
 *
 * The MNA stamps are:
 * 1. KCL at positive node: -I_mic
 * 2. KCL at negative node: +I_mic
 * 3. Constitutive equation (auxiliary row): V_pos - V_neg - R_total_eff * I_mic = V_eff_source
 *
 * @param A The MNA matrix to which the stamps are applied.
 * @param b The MNA right-hand side vector to which the stamps are applied.
 * @param x_current_guess The current guess for node voltages and branch currents (not directly used for this linear model).
 * @param prev_solution The solution from the previous time step, used for the historical term of the inductor companion model.
 * @param time The current simulation time.
 * @param dt The time step size.
 */
void RibbonMicrophone::getStamps(
    Eigen::MatrixXd& A, Eigen::VectorXd& b,
    const Eigen::VectorXd& x_current_guess, const Eigen::VectorXd& prev_solution,
    double time, double dt
) {
    // Get the global (1-based) indices for the microphone's positive and negative nodes.
    // Assuming getNodeIndex(node_name) returns the 1-based index (0 for ground).
    int idx_p_1based = getNodeIndex(pos_node_name);
    int idx_n_1based = getNodeIndex(neg_node_name);

    // Get the global (1-based) index for the auxiliary current variable introduced by this component.
    // Assuming getAuxiliaryVariableStartIndex() returns the 1-based index of the first (and only)
    // auxiliary variable for this microphone component.
    int idx_mic_current_1based = getAuxiliaryVariableStartIndex();

    // --- Calculate terms for the companion model ---

    // Calculate instantaneous sound pressure
    double sound_pressure_at_t = V_sound_pressure_func(time);
    // Calculate voltage generated by the microphone
    double V_mic_generated = sensitivity * sound_pressure_at_t;

    // Calculate effective resistance for the inductor companion model (Backward Euler)
    // For DC analysis (dt=0), an inductor acts as a short circuit, so R_L_eff = 0.
    // For transient analysis (dt > 0), R_L_eff = L_out / dt.
    double R_L_eff = 0.0;
    if (dt > 1e-12) { // Check if dt is actually non-zero for transient analysis
        R_L_eff = L_out / dt;
    }
    // If dt is 0, R_L_eff remains 0, which correctly models the inductor as a short circuit in DC.

    // Get the previous current through the microphone branch from prev_solution.
    // This is the historical term I_mic(t-dt) for the inductor companion model.
    // If prev_solution is not yet populated (e.g., first time step), assume previous current is 0.
    double I_mic_prev = 0.0;
    if (prev_solution.size() > idx_mic_current_1based) {
        I_mic_prev = prev_solution(idx_mic_current_1based);
    }

    // Calculate total effective series resistance of the microphone branch
    // This includes output resistance and effective resistance from the inductor companion model.
    double R_total_eff = R_out + R_L_eff;

    // Calculate the effective voltage source for the companion model.
    // This includes the voltage generated by sound and the historical term from the inductor.
    // V_eff_source = V_mic_generated(t) - (L_out/dt) * I_mic(t-dt)
    double V_L_companion_source = R_L_eff * I_mic_prev; // The voltage source part of the inductor companion model
    double V_eff_source = V_mic_generated - V_L_companion_source;

    // --- Apply MNA stamps to matrix A and vector b ---

    // 1. KCL equation for the positive node (idx_p_1based):
    // Current I_mic leaves the positive node. So, its coefficient is -1.
    // A(row_idx, col_idx) += value
    A(idx_p_1based, idx_mic_current_1based) += -1.0;

    // 2. KCL equation for the negative node (idx_n_1based):
    // Current I_mic enters the negative node. So, its coefficient is +1.
    A(idx_n_1based, idx_mic_current_1based) += 1.0;

    // 3. Constitutive equation for the microphone branch (auxiliary row idx_mic_current_1based):
    // The equation is: V_pos - V_neg - R_total_eff * I_mic = V_eff_source
    // Coefficient for V_pos is +1
    A(idx_mic_current_1based, idx_p_1based) += 1.0;
    // Coefficient for V_neg is -1
    A(idx_mic_current_1based, idx_n_1based) -= 1.0; // Corrected: should be -= 1.0
    // Coefficient for I_mic is -R_total_eff
    A(idx_mic_current_1based, idx_mic_current_1based) -= R_total_eff; // Corrected: should be -= R_total_eff

    // Add the effective voltage source to the right-hand side (RHS) vector b.
    // This term is a known value for the current time step.
    b(idx_mic_current_1based) += V_eff_source;
}

/**
 * @brief Updates the internal state of the RibbonMicrophone.
 *
 * This method updates the previous current through the inductor for the next time step.
 *
 * @param v_curr The current voltage across the component (not directly used for state update here).
 * @param i_curr The current flowing through the component (used to update previous current).
 */
void RibbonMicrophone::updateState(double v_curr, double i_curr) {
    // The previous current for the inductor companion model is the branch current itself.
    // This current is obtained from the MNA solution (x_current_guess) at the auxiliary variable index.
    // However, the Component base class `updateState` is typically called with `v_curr` and `i_curr`
    // which might represent the voltage/current of the *component's main terminals*.
    // For a complex component like this, the actual branch current (I_mic) is the auxiliary variable.
    // So, we need to ensure the solver passes the correct branch current to this `updateState` or
    // access it from the full solution vector `x_current_guess` if that's how it's designed.
    // For now, assuming `i_curr` passed to this function is the branch current I_mic.
    // If not, this needs to be adapted to access x_current_guess(getAuxiliaryVariableStartIndex()).
    // For now, let's assume `i_curr` is the correct branch current.
    // This component does not have an explicit `prev_IL` like in `Cabinet_bass_reflex`,
    // as the historical term is directly handled in `getStamps` using `prev_solution`.
    // So, this `updateState` for RibbonMicrophone might not be strictly necessary
    // if `prev_solution` is always the full solution vector from the previous step.
    // If `prev_solution` in `getStamps` *is* the full solution, then the solver
    // is responsible for passing the correct previous current to `getStamps`.
    // This `updateState` method might be redundant for this specific component
    // if the MNA solver manages `prev_solution` globally.
    // Let's leave it empty as the state update is implicitly handled by the solver
    // passing the previous solution to `getStamps`.
}
