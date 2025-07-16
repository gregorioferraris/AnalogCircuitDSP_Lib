// components/Splitter.h
#ifndef SPLITTER_H
#define SPLITTER_H

#include "Component.h"   // Include the base Component class
#include <Eigen/Dense>   // For Eigen matrix and vector operations
#include <string>        // For std::string

/**
 * @class Splitter
 * @brief Represents a conceptual "splitter" or a very low-resistance connection in a circuit simulation.
 *
 * This class models a component that acts as a near-ideal short circuit between its two nodes.
 * It is implemented as a resistor with a very small, fixed resistance to avoid numerical issues
 * associated with true zero resistance and to fit within the existing MNA framework without
 * introducing auxiliary branch current variables.
 */
class Splitter : public Component {
public:
    // A very small resistance to simulate a near-ideal short circuit.
    // This value is used internally and is not exposed as a constructor parameter
    // to reflect the splitter's conceptual role as a perfect connection.
    const double R_SPLIT = 1e-9; // A very small resistance in Ohms

    /**
     * @brief Constructor for the Splitter component.
     *
     * Initializes the splitter with its name and connected nodes. It represents
     * a conceptual "splitter" or a very low-resistance connection, effectively
     * acting as a near-ideal short circuit. The internal resistance is fixed at R_SPLIT.
     *
     * @param name The unique name of the splitter.
     * @param node1 The name of the first connected node.
     * @param node2 The name of the second connected node.
     * @param tolerance_percent Optional tolerance percentage for component value (not directly used for R_SPLIT).
     */
    Splitter(const std::string& name, const std::string& node1, const std::string& node2,
             double tolerance_percent = 0.0);

    /**
     * @brief Creates a deep copy of the Splitter object.
     * @return A pointer to a new Splitter object, which is a copy of the current instance.
     */
    Component* clone() const override { return new Splitter(*this); }

    /**
     * @brief Applies the stamps of the splitter (as a very low-resistance connection)
     * to the MNA matrix (A) and vector (b).
     *
     * This method treats the splitter as a resistor with a very small, fixed resistance
     * (R_SPLIT) to simulate a near-ideal short circuit. The stamps are applied based
     * on the conductance (G = 1/R_SPLIT).
     *
     * @param A The MNA matrix to which stamps are applied.
     * @param b The MNA right-hand side vector to which stamps are applied.
     * @param x_current_guess The current guess for node voltages and branch currents (not used for this static component).
     * @param prev_solution The solution from the previous time step (not used for this static component).
     * @param time The current simulation time (not used for this static component).
     * @param dt The time step size (not used for this static component).
     */
    void getStamps(
        Eigen::MatrixXd& A, Eigen::VectorXd& b,
        const Eigen::VectorXd& x_current_guess, const Eigen::VectorXd& prev_solution,
        double time, double dt
    ) override;

    /**
     * @brief Updates the internal state of the splitter.
     *
     * As the splitter is modeled as a static, near-ideal connection (its behavior
     * depends only on the instantaneous voltages across it, not on past states),
     * this method does not perform any state updates. It is provided to satisfy
     * the Component base class interface.
     *
     * @param v_curr The current voltage across the component (not used).
     * @param i_curr The current flowing through the component (not used).
     */
    void updateState(double v_curr, double i_curr) override;
};

#endif // SPLITTER_H
