// components/CloseBoxCabinet.h
#ifndef CLOSE_BOX_CABINET_H
#define CLOSE_BOX_CABINET_H

#include "Component.h" // Include the base Component class

/**
 * @brief Represents a Close Box Cabinet component in a circuit simulation.
 *
 * This component typically represents a physical enclosure or a conceptual grouping
 * of other components. By default, it does not introduce any direct electrical
 * stamps to the MNA matrix unless specific parasitic or shielding effects are modeled.
 * It inherits from the base 'Component' class, requiring two nodes for its definition
 * even if they are conceptual for this specific component.
 */
class CloseBoxCabinet : public Component {
public:
    /**
     * @brief Constructor for the CloseBoxCabinet component.
     *
     * Initializes a Close Box Cabinet with its name and two connected nodes.
     *
     * @param name The unique name of the cabinet.
     * @param node1 The name of the first conceptual node associated with the cabinet.
     * @param node2 The name of the second conceptual node associated with the cabinet.
     * @param tolerance_percent Optional tolerance percentage (inherited from Component).
     */
    CloseBoxCabinet(const std::string& name,
                    const std::string& node1, const std::string& node2,
                    double tolerance_percent = 0.0);

    /**
     * @brief Applies the stamps of the CloseBoxCabinet to the MNA matrix (A) and vector (b).
     *
     * For a basic CloseBoxCabinet, this method typically does not add any stamps
     * to the MNA matrix or vector, as it represents a physical enclosure rather
     * than a direct electrical element like a resistor or capacitor.
     * However, it is overridden to comply with the Component interface.
     *
     * If specific electrical effects (e.g., parasitic capacitance to ground,
     * shielding effects) need to be modeled, this method would be extended
     * to include the appropriate stamps.
     *
     * @param A The MNA matrix to which stamps are applied.
     * @param b The MNA right-hand side vector to which stamps are applied.
     * @param x_current_guess The current guess for node voltages and branch currents.
     * @param prev_solution The solution from the previous time step.
     * @param time The current simulation time.
     * @param dt The time step size.
     */
    void getStamps(
        Eigen::MatrixXd& A, Eigen::VectorXd& b,
        const Eigen::VectorXd& x_current_guess, const Eigen::VectorXd& prev_solution,
        double time, double dt
    ) override;

    // No specific state update needed for a passive, non-stamping component
    // If it were to model complex effects, updateState might be needed.
    // void updateState(...) override;
};

#endif // CLOSE_BOX_CABINET_H
