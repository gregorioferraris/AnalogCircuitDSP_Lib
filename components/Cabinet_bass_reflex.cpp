// components/Cabinet_bass_reflex.cpp
#include "Cabinet_bass_reflex.h"
#include <iostream> // For error messages
#include <limits>   // For numeric_limits

/**
 * @brief Constructor for the Cabinet_bass_reflex component.
 *
 * Initializes the bass-reflex cabinet model with its name, connected nodes,
 * tuning frequency, quality factor, and equivalent loading resistance.
 * It also calculates the equivalent L and C values and initializes state variables.
 *
 * @param name The unique name of the component.
 * @param node1 The name of the first connected node.
 * @param node2 The name of the second connected node.
 * @param Fb The tuning frequency of the bass-reflex system (Hz). Must be > 0.
 * @param Qb The quality factor of the bass-reflex system at Fb. Must be > 0.
 * @param R_load The equivalent electrical resistance for damping/loading (Ohms). Must be > 0.
 * @param tolerance_percent Optional tolerance percentage.
 */
Cabinet_bass_reflex::Cabinet_bass_reflex(const std::string& name, const std::string& node1, const std::string& node2,
                                         double Fb, double Qb, double R_load, double tolerance_percent)
    : Component(name, node1, node2, tolerance_percent),
      Fb(Fb), Qb(Qb), R_load(R_load),
      L_eq(0.0), C_eq(0.0), // Will be calculated below
      prev_IL(0.0), prev_VC(0.0) // Initialize state variables for transient analysis
{
    // A parallel RLC equivalent does not typically require additional auxiliary variables
    // beyond the node voltages, as the companion models directly stamp into the nodal equations.
    setNumAuxiliaryVariables(0);

    // Validate input parameters
    if (Fb <= 0.0) {
        std::cerr << "Warning: Cabinet_bass_reflex " << name << " has a non-positive tuning frequency (Fb). Setting to 1.0Hz." << std::endl;
        this->Fb = 1.0;
    }
    if (Qb <= 0.0) {
        std::cerr << "Warning: Cabinet_bass_reflex " << name << " has a non-positive quality factor (Qb). Setting to 0.707." << std::endl;
        this->Qb = 0.707; // Typical value for Butterworth response
    }
    if (R_load <= 0.0) {
        std::cerr << "Warning: Cabinet_bass_reflex " << name << " has a non-positive loading resistance (R_load). Setting to 8.0 Ohms." << std::endl;
        this->R_load = 8.0; // Typical loudspeaker impedance
    }

    // Calculate equivalent L and C for a parallel RLC circuit
    // Angular frequency: omega_b = 2 * PI * Fb
    // For a parallel RLC:
    // Q = R_load * sqrt(C_eq / L_eq)
    // omega_b = 1 / sqrt(L_eq * C_eq)
    //
    // From these, we derive:
    // C_eq = Qb / (R_load * omega_b)
    // L_eq = R_load / (Qb * omega_b)

    double omega_b = 2.0 * M_PI * this->Fb;

    if (omega_b == 0.0) { // Avoid division by zero if Fb was somehow 0 despite check
        std::cerr << "Error: Calculated omega_b is zero for Cabinet_bass_reflex " << name << ". Cannot calculate L_eq/C_eq." << std::endl;
        L_eq = std::numeric_limits<double>::infinity(); // Indicate invalid state
        C_eq = std::numeric_limits<double>::infinity();
    } else {
        C_eq = this->Qb / (this->R_load * omega_b);
        L_eq = this->R_load / (this->Qb * omega_b);
    }

    std::cout << "Cabinet_bass_reflex " << name << " initialized with Fb=" << Fb << "Hz, Qb=" << Qb
              << ", R_load=" << R_load << " Ohms. Equivalent L_eq=" << L_eq << " H, C_eq=" << C_eq << " F." << std::endl;
}

/**
 * @brief Applies the stamps of the equivalent parallel RLC circuit to the MNA matrix (A) and vector (b).
 *
 * This method uses Backward Euler companion models for the inductor and capacitor
 * to stamp their contributions based on the previous time step's state.
 *
 * @param A The MNA matrix to which stamps are applied.
 * @param b The MNA right-hand side vector to which stamps are applied.
 * @param x_current_guess The current guess for node voltages and branch currents (not used for stamping, but for updateState).
 * @param prev_solution The solution from the previous time step (used for companion models).
 * @param time The current simulation time.
 * @param dt The time step size. Crucial for companion models.
 */
void Cabinet_bass_reflex::getStamps(
    Eigen::MatrixXd& A, Eigen::VectorXd& b,
    const Eigen::VectorXd& x_current_guess, const Eigen::VectorXd& prev_solution,
    double time, double dt
) {
    // Get the global indices for the two connected nodes.
    int idx1 = getNodeIndex(node1);
    int idx2 = getNodeIndex(node2);

    // Check for valid indices
    if (idx1 == -1 || idx2 == -1) {
        std::cerr << "Error: Invalid node index for Cabinet_bass_reflex " << name << std::endl;
        return;
    }

    // Ensure dt is positive for transient analysis
    if (dt <= 0.0) {
        std::cerr << "Error: Time step (dt) must be positive for Cabinet_bass_reflex " << name << " transient stamping." << std::endl;
        return;
    }

    // --- Stamp the parallel Resistor (R_load) ---
    // Current through resistor: I_R = (V(idx1) - V(idx2)) / R_load
    double conductance_R = 1.0 / R_load;
    A(idx1, idx1) += conductance_R;
    A(idx2, idx2) += conductance_R;
    A(idx1, idx2) -= conductance_R;
    A(idx2, idx1) -= conductance_R;

    // --- Stamp the parallel Inductor (L_eq) using Backward Euler companion model ---
    // V_L(t) = L_eq * dI_L/dt  =>  V_L(t) = L_eq * (I_L(t) - I_L(t-dt)) / dt
    // Rearranging for I_L(t): I_L(t) = (dt / L_eq) * V_L(t) + I_L(t-dt)
    // This is equivalent to a conductance G_L = dt / L_eq in parallel with a current source I_srcL = I_L(t-dt)
    if (L_eq == 0.0) { // Handle ideal short circuit if L_eq is zero
        // This case implies an infinite frequency response, effectively a short.
        // For practical purposes, a very small L_eq might be used, or this component
        // might not be suitable for such extreme values.
        std::cerr << "Warning: Cabinet_bass_reflex " << name << " has L_eq=0. Treating as short circuit." << std::endl;
        // Effectively add a very large conductance or handle as a short.
        // For now, let's just avoid division by zero and skip inductor stamping.
    } else {
        double conductance_L = dt / L_eq;
        double I_srcL = prev_IL; // Inductor current from previous time step

        A(idx1, idx1) += conductance_L;
        A(idx2, idx2) += conductance_L;
        A(idx1, idx2) -= conductance_L;
        A(idx2, idx1) -= conductance_L;

        b(idx1) -= I_srcL; // Current source flows from node2 to node1
        b(idx2) += I_srcL; // (opposite direction of voltage drop V_L = V(idx1)-V(idx2))
    }

    // --- Stamp the parallel Capacitor (C_eq) using Backward Euler companion model ---
    // I_C(t) = C_eq * dV_C/dt  =>  I_C(t) = C_eq * (V_C(t) - V_C(t-dt)) / dt
    // Rearranging for I_C(t): I_C(t) = (C_eq / dt) * V_C(t) - (C_eq / dt) * V_C(t-dt)
    // This is equivalent to a conductance G_C = C_eq / dt in parallel with a current source I_srcC = (C_eq / dt) * V_C(t-dt)
    if (C_eq == 0.0) { // Handle ideal open circuit if C_eq is zero
        // This effectively means no capacitor current, so no stamp needed.
    } else {
        double conductance_C = C_eq / dt;
        double V_C_prev = prev_VC; // Capacitor voltage from previous time step
        double I_srcC = conductance_C * V_C_prev;

        A(idx1, idx1) += conductance_C;
        A(idx2, idx2) += conductance_C;
        A(idx1, idx2) -= conductance_C;
        A(idx2, idx1) -= conductance_C;

        b(idx1) -= I_srcC; // Current source flows from node2 to node1
        b(idx2) += I_srcC; // (opposite direction of voltage drop V_C = V(idx1)-V(idx2))
    }
}

/**
 * @brief Updates the internal state variables (prev_IL, prev_VC) for the next time step.
 *
 * This method is called after the MNA system is solved for the current time step.
 * It uses the newly calculated voltage across the component (v_curr) and the previous
 * inductor current (prev_IL) to update the state for the next iteration/time step.
 *
 * @param v_curr The current voltage across the component (V(node1) - V(node2)).
 * @param i_curr The current flowing through the component (not directly used for updating internal state, but available).
 */
void Cabinet_bass_reflex::updateState(double v_curr, double i_curr) {
    // Update capacitor voltage for the next time step
    prev_VC = v_curr;

    // Update inductor current for the next time step using Backward Euler formula
    // I_L(t) = I_L(t-dt) + (dt / L_eq) * V_L(t)
    // Here, v_curr is V_L(t)
    if (L_eq != 0.0) { // Avoid division by zero
        prev_IL += (getDt() / L_eq) * v_curr;
    } else {
        // If L_eq is 0, it's a short circuit, current can be infinite or undefined.
        // For practical simulation, this case should be handled by the user providing non-zero L_eq.
        prev_IL = std::numeric_limits<double>::quiet_NaN(); // Indicate invalid state
    }
}
