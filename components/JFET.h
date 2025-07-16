// components/JFET.h
#ifndef JFET_H
#define JFET_H

#include "Component.h" // Include the base Component class
#include <Eigen/Dense> // For Eigen matrix and vector operations
#include <string>      // For std::string
#include <cmath>       // For std::pow

/**
 * @class JFET
 * @brief Represents a Junction Field-Effect Transistor (JFET) using a simplified model.
 *
 * This class models an N-channel JFET as a voltage-controlled current source (VCCS),
 * where the drain current (ID) is controlled by the Gate-Source voltage (VGS).
 *
 * The model uses the common square-law characteristic for the saturation region:
 * ID = Idss * (1 - VGS / Vp)^2
 *
 * For transient analysis and handling non-linearity, a linearized companion model
 * based on Newton-Raphson is employed. This requires calculating the transconductance (gm)
 * and an equivalent current source from the previous time step's VGS.
 *
 * This basic model does not include parasitic capacitances, channel resistances,
 * or gate-source/gate-drain diodes for simplicity.
 */
class JFET : public Component {
public:
    // Specific JFET nodes
    std::string node_gate;  ///< Name of the Gate node.
    std::string node_drain; ///< Name of the Drain node.
    std::string node_source; ///< Name of the Source node.

    // JFET parameters
    double Idss; ///< Drain current with Gate shorted to Source (VGS = 0) (Amps).
    double Vp;   ///< Pinch-off voltage (Volts, negative for N-channel JFETs).

    // State variable for companion model
    double prev_VGS; ///< Gate-Source voltage from the previous time step.

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
     * @param tolerance_percent Optional tolerance percentage (not directly used for this model).
     */
    JFET(const std::string& name,
         const std::string& node_gate, const std::string& node_drain, const std::string& node_source,
         double Idss, double Vp, double tolerance_percent = 0.0);

    /**
     * @brief Creates a deep copy of the JFET object.
     * @return A pointer to a new JFET object, which is a copy of the current instance.
     */
    Component* clone() const override { return new JFET(*this); }

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
    void getStamps(
        Eigen::MatrixXd& A, Eigen::VectorXd& b,
        const Eigen::VectorXd& x_current_guess, const Eigen::VectorXd& prev_solution,
        double time, double dt
    ) override;

    /**
     * @brief Updates the internal state variable (prev_VGS) for the next time step.
     *
     * This method is called after the MNA system is solved for the current time step.
     * It extracts the current Gate-Source voltage from the solution vector and stores it
     * for use in the companion model of the next time step.
     *
     * @param v_curr The current voltage across the component (not directly used for JFET state).
     * @param i_curr The current flowing through the component (not directly used).
     */
    void updateState(double v_curr, double i_curr) override;

private:
    /**
     * @brief Helper function to calculate the drain current (ID) based on VGS.
     *
     * Implements the square-law characteristic: ID = Idss * (1 - VGS / Vp)^2.
     * Handles the cutoff region (VGS >= Vp) where ID = 0.
     *
     * @param VGS The Gate-Source voltage.
     * @return The calculated drain current.
     */
    double calculateDrainCurrent(double VGS);

    /**
     * @brief Helper function to calculate the transconductance (gm) based on VGS.
     *
     * gm = d(ID)/d(VGS) = -2 * Idss / Vp * (1 - VGS / Vp).
     * Handles the cutoff region (VGS >= Vp) where gm = 0.
     *
     * @param VGS The Gate-Source voltage.
     * @return The calculated transconductance.
     */
    double calculateTransconductance(double VGS);
};

#endif // JFET_H
