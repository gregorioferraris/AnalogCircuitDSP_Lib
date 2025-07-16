// components/Led.h
#ifndef LED_H
#define LED_H

#include "Component.h" // Include the base Component class
#include <Eigen/Dense>   // For Eigen matrix and vector operations
#include <string>        // For std::string

/**
 * @class Led
 * @brief Represents a Light Emitting Diode (LED) component in a circuit simulation.
 *
 * This class models an LED using a smoothed piecewise linear approximation for its
 * current-voltage characteristic. It includes parameters for forward voltage (Vf)
 * and series resistance (Rs).
 */
class Led : public Component {
public:
    double Vf;          // Forward voltage (threshold voltage for conduction)
    double Rs;          // Series resistance (resistance in series with the ideal diode)

    /**
     * @brief Constructor for the Led component.
     * @param name The unique name of the LED.
     * @param node1 The name of the anode node.
     * @param node2 The name of the cathode node.
     * @param forward_voltage The nominal forward voltage (Vf) of the LED.
     * @param series_resistance The nominal series resistance (Rs) of the LED.
     * @param tolerance_percent Optional tolerance percentage for component value.
     */
    Led(const std::string& name, const std::string& node1, const std::string& node2,
        double forward_voltage, double series_resistance, double tolerance_percent = 0.0);

    /**
     * @brief Creates a deep copy of the Led object.
     * @return A pointer to a new Led object, which is a copy of the current instance.
     */
    Component* clone() const override { return new Led(*this); }

    /**
     * @brief Applies the stamps of the LED to the MNA (Modified Nodal Analysis) matrix (A) and vector (b).
     *
     * This method calculates the non-linear current and dynamic conductance of the LED
     * based on the current voltage guess and contributes them to the MNA system.
     * It uses a smoothed exponential model to approximate the LED's I-V characteristic.
     *
     * @param A The MNA matrix to which stamps are applied.
     * @param b The MNA right-hand side vector to which stamps are applied.
     * @param x_current_guess The current guess for node voltages and branch currents (for Newton-Raphson).
     * @param prev_solution The solution from the previous time step (not directly used for static LED model).
     * @param time The current simulation time (not directly used for static LED model).
     * @param dt The time step size (not directly used for static LED model).
     */
    void getStamps(
        Eigen::MatrixXd& A, Eigen::VectorXd& b,
        const Eigen::VectorXd& x_current_guess, const Eigen::VectorXd& prev_solution,
        double time, double dt
    ) override;

    /**
     * @brief Updates the internal state of the LED.
     *
     * For a static non-linear component like an LED, there is typically no internal
     * state that evolves over time based on its own voltage/current history.
     * This method is left empty for this model.
     *
     * @param v_curr The current voltage across the component.
     * @param i_curr The current flowing through the component.
     */
    void updateState(double v_curr, double i_curr) override;
};

#endif // LED_H
