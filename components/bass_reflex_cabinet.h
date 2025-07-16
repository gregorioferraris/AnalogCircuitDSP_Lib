// components/Cabinet_bass_reflex.h
#ifndef CABINET_BASS_REFLEX_H
#define CABINET_BASS_REFLEX_H

#include "Component.h" // Include the base Component class
#include <Eigen/Dense> // For Eigen matrix and vector operations
#include <string>      // For std::string
#include <cmath>       // For M_PI (pi constant)

/**
 * @class Cabinet_bass_reflex
 * @brief Represents an equivalent electrical model of a bass-reflex loudspeaker cabinet.
 *
 * This class models the acoustic impedance of a bass-reflex cabinet and port
 * as an equivalent parallel RLC circuit in the electrical domain. This allows
 * its integration into a circuit simulation using Modified Nodal Analysis (MNA).
 *
 * The component connects between two electrical nodes (node1 and node2) and
 * applies stamps for its equivalent RLC circuit. For transient analysis,
 * companion models (Backward Euler) are used for the inductor and capacitor,
 * requiring storage of the previous time-step's inductor current and capacitor voltage.
 *
 * The RLC parameters (L_eq, C_eq, R_load) are derived from the specified
 * tuning frequency (Fb), quality factor (Qb), and an equivalent electrical
 * loading resistance (R_load).
 *
 * The equivalent parallel RLC circuit parameters are calculated as:
 * Angular frequency: omega_b = 2 * PI * Fb
 * Equivalent Capacitance: C_eq = Qb / (R_load * omega_b)
 * Equivalent Inductance: L_eq = R_load / (Qb * omega_b)
 */
class Cabinet_bass_reflex : public Component {
public:
    double Fb;       ///< Tuning frequency of the bass-reflex system (Hz).
    double Qb;       ///< Quality factor of the bass-reflex system at Fb.
    double R_load;   ///< Equivalent electrical resistance representing the damping/loading (Ohms).

    double L_eq;     ///< Calculated equivalent inductance (Henries).
    double C_eq;     ///< Calculated equivalent capacitance (Farads).

    // State variables for companion models
    double prev_IL;  ///< Inductor current from the previous time step.
    double prev_VC;  ///< Capacitor voltage from the previous time step.

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
     * @param tolerance_percent Optional tolerance percentage for component value (not directly used for this complex model).
     */
    Cabinet_bass_reflex(const std::string& name, const std::string& node1, const std::string& node2,
                        double Fb, double Qb, double R_load, double tolerance_percent = 0.0);

    /**
     * @brief Creates a deep copy of the Cabinet_bass_reflex object.
     * @return A pointer to a new Cabinet_bass_reflex object, which is a copy of the current instance.
     */
    Component* clone() const override { return new Cabinet_bass_reflex(*this); }

    /**
     * @brief Applies the stamps of the equivalent parallel RLC circuit to the MNA matrix (A) and vector (b).
     *
     * This method uses Backward Euler companion models for the inductor and capacitor
     * to stamp their contributions based on the previous time step's state.
     *
     * @param A The MNA matrix to which stamps are applied.
     * @param b The MNA right-hand side vector to which stamps are applied.
     * @param x_current_guess The current guess for node voltages and branch currents (used for updating state).
     * @param prev_solution The solution from the previous time step (used for companion models).
     * @param time The current simulation time (not directly used for stamping).
     * @param dt The time step size. Crucial for companion models.
     */
    void getStamps(
        Eigen::MatrixXd& A, Eigen::VectorXd& b,
        const Eigen::VectorXd& x_current_guess, const Eigen::VectorXd& prev_solution,
        double time, double dt
    ) override;

    /**
     * @brief Updates the internal state variables (prev_IL, prev_VC) for the next time step.
     *
     * This method is called after the MNA system is solved for the current time step.
     * It uses the newly calculated voltage across the component and the previous
     * inductor current to update the state for the next iteration/time step.
     *
     * @param v_curr The current voltage across the component (V(node1) - V(node2)).
     * @param i_curr The current flowing through the component (not directly used for updating internal state, but available).
     */
    void updateState(double v_curr, double i_curr) override;
};

#endif // CABINET_BASS_REFLEX_H
