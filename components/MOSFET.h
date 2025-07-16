// components/MOSFET.h
#ifndef MOSFET_H
#define MOSFET_H

#include "Component.h" // Include the base Component class
#include <string>
#include <vector>
#include <cmath> // For std::sqrt, std::pow

/**
 * @brief Represents a MOSFET (Metal-Oxide-Semiconductor Field-Effect Transistor)
 * component in a circuit simulation, using a simplified SPICE Level 1 model.
 *
 * This class models an N-channel MOSFET (NMOS). For PMOS, parameters and
 * voltage polarities would need to be adjusted or a separate PMOS class created.
 *
 * It has four terminals: Drain (D), Gate (G), Source (S), and Bulk/Body (B).
 * The companion model for non-linear devices is implemented in getStamps.
 */
class MOSFET : public Component {
public:
    /**
     * @brief Constructor for the MOSFET component.
     *
     * Initializes a MOSFET with its name, four connected nodes (Drain, Gate, Source, Bulk),
     * and key SPICE Level 1 model parameters.
     *
     * @param name The unique name of the MOSFET.
     * @param drain_node The name of the Drain node.
     * @param gate_node The name of the Gate node.
     * @param source_node The name of the Source node.
     * @param bulk_node The name of the Bulk/Body node.
     * @param type The type of MOSFET ("NMOS" or "PMOS").
     * @param kp Transconductance parameter (k' = mu_n * C_ox) in A/V^2.
     * @param vt0 Zero-bias threshold voltage (Vth0) in V.
     * @param lambda Channel-length modulation parameter (lambda) in V^-1.
     * @param gamma Bulk threshold parameter (gamma) in V^(1/2).
     * @param phi Surface potential (phi) in V.
     * @param W Channel width in meters.
     * @param L Channel length in meters.
     * @param tolerance_percent Optional tolerance percentage (inherited from Component).
     */
    MOSFET(const std::string& name,
           const std::string& drain_node, const std::string& gate_node,
           const std::string& source_node, const std::string& bulk_node,
           const std::string& type,
           double kp, double vt0, double lambda, double gamma, double phi,
           double W, double L,
           double tolerance_percent = 0.0);

    /**
     * @brief Applies the stamps of the MOSFET companion model to the MNA matrix (A) and vector (b).
     *
     * This method implements the non-linear companion model for the MOSFET.
     * It calculates the drain current (Id) and small-signal conductances (gds, gm, gmb)
     * based on the current guess of node voltages (x_current_guess).
     * These values are then used to stamp the MNA matrix and vector for the current iteration.
     *
     * @param A The MNA matrix to which stamps are applied.
     * @param b The MNA right-hand side vector to which stamps are applied.
     * @param x_current_guess The current guess for node voltages and branch currents.
     * @param prev_solution The solution from the previous time step (not directly used in this non-linear model for time-stepping, but could be for dynamic models).
     * @param time The current simulation time.
     * @param dt The time step size.
     */
    void getStamps(
        Eigen::MatrixXd& A, Eigen::VectorXd& b,
        const Eigen::VectorXd& x_current_guess, const Eigen::VectorXd& prev_solution,
        double time, double dt
    ) override;

private:
    // MOSFET model parameters
    std::string type; // "NMOS" or "PMOS"
    double kp;        // Transconductance parameter (k' = mu_n * C_ox)
    double vt0;       // Zero-bias threshold voltage
    double lambda;    // Channel-length modulation parameter
    double gamma;     // Bulk threshold parameter
    double phi;       // Surface potential (2*phi_F)
    double W;         // Channel width
    double L;         // Channel length
    double Leff;      // Effective channel length (L - 2*LD, if LD is a parameter)

    // Node names for convenience
    std::string drain_node_name;
    std::string gate_node_name;
    std::string source_node_name;
    std::string bulk_node_name;

    /**
     * @brief Calculates the threshold voltage (Vth) considering the body effect.
     * @param Vsb The source-bulk voltage.
     * @return The threshold voltage.
     */
    double calculateThresholdVoltage(double Vsb) const;

    /**
     * @brief Calculates the drain current (Id) based on operating region.
     * @param Vgs Gate-Source voltage.
     * @param Vds Drain-Source voltage.
     * @param Vth Threshold voltage.
     * @return The drain current.
     */
    double calculateDrainCurrent(double Vgs, double Vds, double Vth) const;

    /**
     * @brief Calculates the output conductance (gds = d(Id)/d(Vds)).
     * @param Vgs Gate-Source voltage.
     * @param Vds Drain-Source voltage.
     * @param Vth Threshold voltage.
     * @return The output conductance.
     */
    double calculateGDS(double Vgs, double Vds, double Vth) const;

    /**
     * @brief Calculates the transconductance (gm = d(Id)/d(Vgs)).
     * @param Vgs Gate-Source voltage.
     * @param Vds Drain-Source voltage.
     * @param Vth Threshold voltage.
     * @return The transconductance.
     */
    double calculateGM(double Vgs, double Vds, double Vth) const;

    /**
     * @brief Calculates the bulk transconductance (gmb = d(Id)/d(Vbs)).
     * This requires the derivative of Vth with respect to Vbs.
     * @param Vgs Gate-Source voltage.
     * @param Vds Drain-Source voltage.
     * @param Vbs Bulk-Source voltage.
     * @param Vth Threshold voltage.
     * @return The bulk transconductance.
     */
    double calculateGMB(double Vgs, double Vds, double Vbs, double Vth) const;
};

#endif // MOSFET_H
