// components/Transformer.h
#ifndef TRANSFORMER_H
#define TRANSFORMER_H

#include "Component.h" // Include the base Component class
#include <Eigen/Dense> // For Eigen matrix and vector operations
#include <string>      // For std::string
#include <vector>      // For std::vector (to store multiple node names if needed)

/**
 * @class Transformer
 * @brief Represents an ideal transformer in a circuit simulation.
 *
 * This class models an ideal transformer, a four-terminal device with a primary
 * and a secondary winding. Its behavior is defined by a turns ratio 'a' (Np/Ns),
 * which relates the voltages and currents between the primary and secondary sides.
 *
 * It introduces two auxiliary variables into the MNA matrix: one for the primary
 * current (Ip) and one for the secondary current (Is).
 *
 * The equations implemented are:
 * 1. V_primary - a * V_secondary = 0
 * 2. Ip + (1/a) * Is = 0
 */
class Transformer : public Component {
public:
    // Node names for the secondary winding
    std::string node1_s; ///< Name of the first node on the secondary side (positive terminal).
    std::string node2_s; ///< Name of the second node on the secondary side (negative terminal).

    double turns_ratio; ///< The turns ratio (Np/Ns) of the transformer.

    /**
     * @brief Constructor for the Transformer component.
     *
     * Initializes an ideal transformer with its name, four connected nodes,
     * and its turns ratio.
     *
     * @param name The unique name of the transformer.
     * @param node1_p The name of the first node on the primary side (positive terminal).
     * @param node2_p The name of the second node on the primary side (negative terminal).
     * @param node1_s The name of the first node on the secondary side (positive terminal).
     * @param node2_s The name of the second node on the secondary side (negative terminal).
     * @param turns_ratio The turns ratio (Np/Ns) of the transformer. Must be non-zero.
     * @param tolerance_percent Optional tolerance percentage for component value (not directly used for ideal transformer).
     */
    Transformer(const std::string& name, const std::string& node1_p, const std::string& node2_p,
                const std::string& node1_s, const std::string& node2_s,
                double turns_ratio, double tolerance_percent = 0.0);

    /**
     * @brief Creates a deep copy of the Transformer object.
     * @return A pointer to a new Transformer object, which is a copy of the current instance.
     */
    Component* clone() const override { return new Transformer(*this); }

    /**
     * @brief Applies the stamps of the ideal transformer to the MNA matrix (A) and vector (b).
     *
     * This method introduces two auxiliary variables for the primary (Ip) and secondary (Is)
     * currents. It stamps the voltage and current relationships of an ideal transformer.
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
     * @brief Updates the internal state of the transformer.
     *
     * As an ideal transformer is a static, lossless component (its behavior
     * depends only on instantaneous voltages and currents, not on past states),
     * this method does not perform any state updates. It is provided to satisfy
     * the Component base class interface.
     *
     * @param v_curr The current voltage across the component (not used).
     * @param i_curr The current flowing through the component (not used).
     */
    void updateState(double v_curr, double i_curr) override;
};

#endif // TRANSFORMER_H
