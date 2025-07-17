// components/Led.cpp
#include "Led.h"
#include <cmath>   // For std::exp and std::log
#include <iostream> // For std::cerr (warning messages)

// Define constants for the LED model's numerical behavior
// LED_SHARPNESS_FACTOR (k): Controls how sharply the LED turns on.
// A higher value makes the transition from off to on more abrupt.
const double LED_SHARPNESS_FACTOR = 500.0;

// LED_LEAKAGE_CURRENT: A small current added for numerical stability,
// especially when the LED is reverse-biased or off. Prevents issues
// with log(0) or extremely small numbers in calculations.
const double LED_LEAKAGE_CURRENT = 1e-12; // Amperes

/**
 * @brief Constructor for the Led component.
 *
 * Initializes the LED with its name, connected nodes, forward voltage,
 * and series resistance. It also performs a basic check to ensure
 * the series resistance is not zero, setting a minimum if necessary.
 *
 * @param name The unique name of the LED.
 * @param node1 The name of the anode node.
 * @param node2 The name of the cathode node.
 * @param forward_voltage The nominal forward voltage (Vf) of the LED.
 * @param series_resistance The nominal series resistance (Rs) of the LED.
 * @param tolerance_percent Optional tolerance percentage for component value.
 */
Led::Led(const std::string& name, const std::string& node1, const std::string& node2,
         double forward_voltage, double series_resistance, double tolerance_percent)
    : Component(name, node1, node2, tolerance_percent), // Call base class constructor
      Vf(forward_voltage),
      Rs(series_resistance)
{
    // Ensure series resistance is not zero or negative to prevent division by zero
    // or unphysical behavior in the model.
    if (Rs <= 0) {
        Rs = 1e-9; // Set a small, non-zero resistance
        std::cerr << "Warning: LED series resistance for component '" << name
                  << "' cannot be zero or negative. Setting to " << Rs << " Ohm." << std::endl;
    }
}

/**
 * @brief Applies the stamps of the LED to the MNA matrix (A) and vector (b).
 *
 * This method implements a smoothed piecewise linear model for the LED's
 * current-voltage characteristic, suitable for Newton-Raphson iteration.
 * The current (I_LED) and dynamic conductance (g_LED) are calculated based
 * on the current voltage guess across the LED.
 *
 * The model used is:
 * I_LED = (1/Rs) * log(1 + exp(k * (V_LED_guess - Vf))) + I_leakage
 * where:
 * k = LED_SHARPNESS_FACTOR
 * V_LED_guess = current voltage across the LED (anode - cathode)
 * Vf = forward voltage
 * Rs = series resistance
 * I_leakage = LED_LEAKAGE_CURRENT
 *
 * The dynamic conductance g_LED is the derivative of I_LED with respect to V_LED_guess:
 * g_LED = (k / Rs) * (exp(k * (V_LED_guess - Vf))) / (1 + exp(k * (V_LED_guess - Vf)))
 *
 * These values are then used to stamp the A matrix and b vector for the MNA solver.
 *
 * @param A The MNA matrix to which stamps are applied.
 * @param b The MNA right-hand side vector to which stamps are applied.
 * @param x_current_guess The current guess for node voltages and branch currents.
 * @param prev_solution The solution from the previous time step (not used by this static model).
 * @param time The current simulation time (not used by this static model).
 * @param dt The time step size (not used by this static model).
 */
void Led::getStamps(
    Eigen::MatrixXd& A, Eigen::VectorXd& b,
    const Eigen::VectorXd& x_current_guess, const Eigen::VectorXd& prev_solution,
    double time, double dt
) {
    // Get the global indices for the anode and cathode nodes.
    // getNodeIndex is assumed to be a method of the base Component class
    // that correctly maps node names to their corresponding indices in the MNA matrix.
    int idx1 = getNodeIndex(node1); // Anode node index
    int idx2 = getNodeIndex(node2); // Cathode node index

    // Handle ground node (index 0) which has voltage 0
    double V1_curr = (idx1 == 0) ? 0.0 : x_current_guess(idx1 - 1);
    double V2_curr = (idx2 == 0) ? 0.0 : x_current_guess(idx2 - 1);

    // Calculate the voltage across the LED (V_anode - V_cathode) based on the current guess.
    double V_LED_guess = V1_curr - V2_curr;

    // Calculate the argument for the exponential function in the smoothed model.
    double arg = LED_SHARPNESS_FACTOR * (V_LED_guess - Vf);

    // Calculate the log term for the current equation using a numerically stable approach.
    // std::log(1.0 + std::exp(x)) can be unstable for very large x.
    // If x > 0, log(1 + exp(x)) = x + log(1 + exp(-x))
    // If x <= 0, log(1 + exp(x))
    double log_term;
    if (arg > 0) {
        log_term = arg + std::log(1.0 + std::exp(-arg));
    } else {
        log_term = std::log(1.0 + std::exp(arg));
    }

    // Calculate the current flowing through the LED.
    // This is the non-linear current I_NL(V_LED_guess).
    double I_LED = (1.0 / Rs) * log_term + LED_LEAKAGE_CURRENT;

    // Calculate the dynamic conductance (g_LED), which is the derivative of I_LED
    // with respect to V_LED_guess. This is used for the Jacobian matrix (A).
    // The derivative of log(1 + exp(x)) is exp(x) / (1 + exp(x)), which is the logistic sigmoid function.
    double exp_arg = std::exp(arg); // Compute exp(arg) once to avoid redundant calculations
    double g_LED = (LED_SHARPNESS_FACTOR / Rs) * (exp_arg / (1.0 + exp_arg));

    // Apply the stamps to the MNA matrix A and vector b.
    // These stamps represent the linearized contribution of the non-linear LED.
    // The general form for a non-linear current I_NL(V_1 - V_2) is:
    // KCL at node 1: +I_NL
    // KCL at node 2: -I_NL
    // Companion model: I_NL = g_NL * (V1 - V2) + I_eq
    // where I_eq = I_NL(V_guess) - g_NL * (V1_guess - V2_guess)

    // KCL at anode node (idx1)
    if (idx1 != 0) {
        A(idx1, idx1) += g_LED;
        if (idx2 != 0) A(idx1, idx2) -= g_LED;
        b(idx1) -= (I_LED - g_LED * V_LED_guess);
    }

    // KCL at cathode node (idx2)
    if (idx2 != 0) {
        A(idx2, idx2) += g_LED;
        if (idx1 != 0) A(idx2, idx1) -= g_LED;
        b(idx2) += (I_LED - g_LED * V_LED_guess);
    }
}

/**
 * @brief Updates the internal state of the LED.
 *
 * For this static non-linear model, there is no internal state to update
 * based on previous time steps. This method is empty.
 *
 * @param v_curr The current voltage across the component (not used).
 * @param i_curr The current flowing through the component (not used).
 */
void Led::updateState(double v_curr, double i_curr) {
    // No state to update for this static non-linear model.
}
