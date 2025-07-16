// components/MOSFET.cpp
#include "MOSFET.h"
#include <iostream> // For error messages and debugging output
#include <algorithm> // For std::max, std::min

// Define a small epsilon for numerical stability, especially for derivatives near boundaries
const double EPSILON = 1e-9; // Small voltage difference for numerical stability

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
MOSFET::MOSFET(const std::string& name,
               const std::string& drain_node, const std::string& gate_node,
               const std::string& source_node, const std::string& bulk_node,
               const std::string& type,
               double kp, double vt0, double lambda, double gamma, double phi,
               double W, double L,
               double tolerance_percent)
    : Component(name, drain_node, gate_node, tolerance_percent), // Using drain and gate as primary nodes for Component base
      type(type), kp(kp), vt0(vt0), lambda(lambda), gamma(gamma), phi(phi), W(W), L(L),
      drain_node_name(drain_node), gate_node_name(gate_node),
      source_node_name(source_node), bulk_node_name(bulk_node)
{
    // MOSFETs do not typically add auxiliary variables for their current,
    // as their non-linear current source is directly integrated into KCL equations.
    setNumAuxiliaryVariables(0);

    // Effective channel length (assuming no lateral diffusion LD for simplicity in Level 1)
    // If LD was a parameter, Leff = L - 2 * LD;
    Leff = L; // For simplicity, assume Leff = L if LD is not provided.

    // Basic validation for parameters
    if (kp <= 0 || W <= 0 || L <= 0 || phi <= 0) {
        std::cerr << "Warning: MOSFET " << name << " has invalid (non-positive) model parameters (kp, W, L, phi)." << std::endl;
        // Consider setting default valid values or throwing an exception
    }
    if (type != "NMOS" && type != "PMOS") {
        std::cerr << "Warning: MOSFET " << name << " has an invalid type. Must be 'NMOS' or 'PMOS'." << std::endl;
    }

    std::cout << "MOSFET " << name << " initialized (" << type << ", W/L=" << W/L << "). D:" << drain_node << " G:" << gate_node << " S:" << source_node << " B:" << bulk_node << std::endl;
}

/**
 * @brief Calculates the threshold voltage (Vth) considering the body effect.
 * Vth = Vth0 + gamma * (sqrt(2*phi + Vsb) - sqrt(2*phi))
 * @param Vsb The source-bulk voltage.
 * @return The threshold voltage.
 */
double MOSFET::calculateThresholdVoltage(double Vsb) const {
    if (type == "NMOS") {
        // For NMOS, Vsb should be >= 0 for sqrt to be real. If Vsb < 0, it's forward biased,
        // but the model typically assumes Vsb >= 0.
        // Clamp Vsb to prevent sqrt of negative numbers, though ideally Vsb should be positive for NMOS.
        double Vsb_clamped = std::max(0.0, Vsb);
        return vt0 + gamma * (std::sqrt(2 * phi + Vsb_clamped) - std::sqrt(2 * phi));
    } else { // PMOS
        // For PMOS, Vsb is typically negative. Vbs is positive.
        // Vth = Vth0 - gamma * (sqrt(2*phi - Vbs) - sqrt(2*phi))
        // Or, using Vsb: Vth = Vth0 + gamma * (sqrt(2*phi - Vsb) - sqrt(2*phi))
        // The body effect for PMOS is usually Vth = Vth0 - gamma * (sqrt(2*phi - Vbs) - sqrt(2*phi))
        // where Vbs is positive. If Vsb is used, Vsb is negative.
        // Using abs(Vbs) or -Vsb for PMOS body effect.
        double Vbs_abs = std::abs(Vsb); // Vbs = -Vsb
        double Vbs_clamped = std::max(0.0, Vbs_abs); // Ensure positive for sqrt
        return vt0 - gamma * (std::sqrt(2 * phi + Vbs_clamped) - std::sqrt(2 * phi));
    }
}

/**
 * @brief Calculates the drain current (Id) based on operating region (SPICE Level 1).
 * Assumes NMOS for calculations, PMOS would have inverted polarities.
 * @param Vgs Gate-Source voltage.
 * @param Vds Drain-Source voltage.
 * @param Vth Threshold voltage.
 * @return The drain current.
 */
double MOSFET::calculateDrainCurrent(double Vgs, double Vds, double Vth) const {
    double id = 0.0;
    double Vov = Vgs - Vth; // Overdrive voltage

    if (type == "NMOS") {
        if (Vgs <= Vth + EPSILON) { // Cutoff region
            id = 0.0;
        } else if (Vds < Vov - EPSILON) { // Triode (Linear) region
            id = kp * (W / Leff) * (Vov * Vds - 0.5 * Vds * Vds) * (1.0 + lambda * Vds);
        } else { // Saturation region (Vds >= Vov)
            id = 0.5 * kp * (W / Leff) * Vov * Vov * (1.0 + lambda * Vds);
        }
    } else { // PMOS
        // For PMOS, Vgs, Vds, Vth are typically negative or compared with negative values.
        // It's often easier to work with absolute values or invert the logic.
        // Let's use |Vgs|, |Vds|, |Vth| and invert the current direction.
        double abs_Vgs = std::abs(Vgs);
        double abs_Vds = std::abs(Vds);
        double abs_Vth = std::abs(Vth); // Vth is usually negative for PMOS

        // PMOS is ON when Vgs < Vth (i.e., |Vgs| > |Vth|)
        // PMOS Triode when |Vds| < |Vgs - Vth|
        // PMOS Saturation when |Vds| >= |Vgs - Vth|

        if (Vgs >= Vth - EPSILON) { // Cutoff region (Vgs less negative than Vth)
            id = 0.0;
        } else if (Vds > Vov + EPSILON) { // Triode (Linear) region (Vds less negative than Vds_sat)
            // Vds is negative, Vov is negative.
            // Id = kp * (W/Leff) * (Vov * Vds - 0.5 * Vds * Vds) * (1.0 - lambda * Vds)
            // Using abs values: Id = 0.5 * kp * (W/Leff) * (2*abs(Vov)*abs(Vds) - abs(Vds)*abs(Vds)) * (1.0 + lambda * abs(Vds))
            // The SPICE model for PMOS usually uses the same equations but with |Vgs|, |Vds|, |Vth|
            // and then multiplies by -1 for current direction.
            id = kp * (W / Leff) * (Vov * Vds - 0.5 * Vds * Vds) * (1.0 - lambda * Vds); // Vds is negative, lambda*Vds becomes negative
        } else { // Saturation region (Vds <= Vov)
            id = 0.5 * kp * (W / Leff) * Vov * Vov * (1.0 - lambda * Vds);
        }
    }
    return id;
}

/**
 * @brief Calculates the output conductance (gds = d(Id)/d(Vds)).
 * @param Vgs Gate-Source voltage.
 * @param Vds Drain-Source voltage.
 * @param Vth Threshold voltage.
 * @return The output conductance.
 */
double MOSFET::calculateGDS(double Vgs, double Vds, double Vth) const {
    double gds = 0.0;
    double Vov = Vgs - Vth;

    if (type == "NMOS") {
        if (Vgs <= Vth + EPSILON) { // Cutoff
            gds = 0.0;
        } else if (Vds < Vov - EPSILON) { // Triode
            gds = kp * (W / Leff) * (Vov - Vds + lambda * (2 * Vov * Vds - 1.5 * Vds * Vds));
        } else { // Saturation
            gds = 0.5 * kp * (W / Leff) * Vov * Vov * lambda;
        }
    } else { // PMOS
        if (Vgs >= Vth - EPSILON) { // Cutoff
            gds = 0.0;
        } else if (Vds > Vov + EPSILON) { // Triode
            // Derivative of kp * (W/Leff) * (Vov * Vds - 0.5 * Vds * Vds) * (1.0 - lambda * Vds) w.r.t Vds
            gds = kp * (W / Leff) * (Vov - Vds - lambda * (Vov * Vds - 0.5 * Vds * Vds) - (1.0 - lambda * Vds) * Vds);
            gds = kp * (W / Leff) * (Vov - Vds + lambda * (Vov * Vds - 0.5 * Vds * Vds) - lambda * (Vov - Vds)); // Simplified
            gds = kp * (W / Leff) * (Vov - Vds - lambda * Vov * Vds + 0.5 * lambda * Vds * Vds + lambda * Vds * Vds - lambda * Vov * Vds); // More detail
            gds = kp * (W / Leff) * (Vov - Vds + lambda * (0.5 * Vds * Vds - Vov * Vds)); // Corrected for PMOS Vds negative
            gds = kp * (W / Leff) * (Vov - Vds + lambda * (2 * Vov * Vds - 1.5 * Vds * Vds)); // Using absolute values for consistency
            gds = kp * (W / Leff) * (Vov - Vds + lambda * (2 * Vov * Vds - 1.5 * Vds * Vds)); // This derivative should be for PMOS
            // Let's re-derive for PMOS: Id = K * (Vov * Vds - 0.5 * Vds^2) * (1 - L * Vds)
            // dId/dVds = K * [ (Vov - Vds)(1 - L*Vds) + (Vov*Vds - 0.5*Vds^2)(-L) ]
            // = K * [ Vov - L*Vov*Vds - Vds + L*Vds^2 - L*Vov*Vds + 0.5*L*Vds^2 ]
            // = K * [ Vov - Vds - 2*L*Vov*Vds + 1.5*L*Vds^2 ]
            gds = kp * (W / Leff) * (Vov - Vds - 2 * lambda * Vov * Vds + 1.5 * lambda * Vds * Vds);
        } else { // Saturation
            gds = 0.5 * kp * (W / Leff) * Vov * Vov * (-lambda); // Note: -lambda for PMOS
        }
    }
    return gds;
}

/**
 * @brief Calculates the transconductance (gm = d(Id)/d(Vgs)).
 * @param Vgs Gate-Source voltage.
 * @param Vds Drain-Source voltage.
 * @param Vth Threshold voltage.
 * @return The transconductance.
 */
double MOSFET::calculateGM(double Vgs, double Vds, double Vth) const {
    double gm = 0.0;
    double Vov = Vgs - Vth;

    if (type == "NMOS") {
        if (Vgs <= Vth + EPSILON) { // Cutoff
            gm = 0.0;
        } else if (Vds < Vov - EPSILON) { // Triode
            gm = kp * (W / Leff) * Vds * (1.0 + lambda * Vds);
        } else { // Saturation
            gm = kp * (W / Leff) * Vov * (1.0 + lambda * Vds);
        }
    } else { // PMOS
        if (Vgs >= Vth - EPSILON) { // Cutoff
            gm = 0.0;
        } else if (Vds > Vov + EPSILON) { // Triode
            // Derivative of kp * (W/Leff) * (Vov * Vds - 0.5 * Vds * Vds) * (1.0 - lambda * Vds) w.r.t Vgs
            // Only Vov depends on Vgs. dVov/dVgs = 1.
            gm = kp * (W / Leff) * Vds * (1.0 - lambda * Vds);
        } else { // Saturation
            gm = kp * (W / Leff) * Vov * (1.0 - lambda * Vds);
        }
    }
    return gm;
}

/**
 * @brief Calculates the bulk transconductance (gmb = d(Id)/d(Vbs)).
 * This requires the derivative of Vth with respect to Vbs.
 * gmb = -gm * d(Vth)/d(Vsb)
 * @param Vgs Gate-Source voltage.
 * @param Vds Drain-Source voltage.
 * @param Vbs Bulk-Source voltage.
 * @param Vth Threshold voltage.
 * @return The bulk transconductance.
 */
double MOSFET::calculateGMB(double Vgs, double Vds, double Vbs, double Vth) const {
    double gmb = 0.0;
    double gm_val = calculateGM(Vgs, Vds, Vth); // Get current gm value

    // Derivative of Vth with respect to Vsb (dVth/dVsb)
    // dVth/dVsb = gamma / (2 * sqrt(2*phi + Vsb))
    // For NMOS, Vsb = -Vbs. So dVth/dVbs = -dVth/dVsb.
    // For PMOS, Vsb is negative, Vbs is positive. Vsb = -Vbs.
    // The body effect term is gamma * (sqrt(2*phi + abs(Vsb)) - sqrt(2*phi))
    // So d(sqrt(2*phi + abs(Vsb)))/d(Vbs) = d(sqrt(2*phi + Vbs))/d(Vbs) = 1 / (2 * sqrt(2*phi + Vbs))
    // For NMOS: Vsb is typically positive. Vbs is negative.
    // dVth/dVbs = -gamma / (2 * sqrt(2*phi + Vsb))
    // For PMOS: Vsb is typically negative. Vbs is positive.
    // dVth/dVbs = gamma / (2 * sqrt(2*phi + Vbs))
    // Let's use Vsb directly for the derivative, and then relate to Vbs.
    // Vsb is Source-Bulk voltage. Vbs is Bulk-Source voltage. Vsb = -Vbs.

    double dVth_dVsb_magnitude = 0.0;
    if (std::abs(2 * phi + Vsb) > EPSILON) { // Avoid division by zero
        dVth_dVsb_magnitude = gamma / (2.0 * std::sqrt(std::abs(2 * phi + Vsb)));
    }

    if (type == "NMOS") {
        // For NMOS, Vsb is usually positive. If Vsb < 0 (forward bias), model might break.
        // gmb = gm * dVth/dVsb * (-1) (since Vbs = -Vsb)
        // dVth/dVsb = gamma / (2 * sqrt(2*phi + Vsb))
        // gmb = gm * gamma / (2 * sqrt(2*phi + Vsb))
        gmb = gm_val * dVth_dVsb_magnitude; // gmb = gm * dVth/dVbs_magnitude
    } else { // PMOS
        // For PMOS, Vsb is usually negative. Vbs is positive.
        // Vth = Vth0 - gamma * (sqrt(2*phi + Vbs) - sqrt(2*phi))
        // dVth/dVbs = -gamma / (2 * sqrt(2*phi + Vbs))
        // gmb = -gm * dVth/dVbs = -gm * (-gamma / (2 * sqrt(2*phi + Vbs))) = gm * gamma / (2 * sqrt(2*phi + Vbs))
        gmb = gm_val * dVth_dVsb_magnitude; // gmb = gm * dVth/dVbs_magnitude
    }

    return gmb;
}


/**
 * @brief Applies the stamps of the MOSFET companion model to the MNA matrix (A) and vector (b).
 *
 * This method implements the non-linear companion model for the MOSFET.
 * It calculates the drain current (Id) and small-signal conductances (gds, gm, gmb)
 * based on the current guess of node voltages (x_current_guess).
 * These values are then used to stamp the MNA matrix and vector for the current iteration.
 *
 * The companion model for a non-linear current source $I_D(V_{GS}, V_{DS}, V_{BS})$ is:
 * $I_D \approx I_D(V_{GS0}, V_{DS0}, V_{BS0}) + g_m \Delta V_{GS} + g_{ds} \Delta V_{DS} + g_{mb} \Delta V_{BS}$
 * where $\Delta V = V - V_0$.
 * Rearranging for MNA:
 * $I_D - g_m V_{GS} - g_{ds} V_{DS} - g_{mb} V_{BS} = I_D(V_0) - g_m V_{GS0} - g_{ds} V_{DS0} - g_{mb} V_{BS0}$
 *
 * So, the terms to stamp are:
 * - Conductances: $g_m$, $g_{ds}$, $g_{mb}$
 * - Equivalent current source: $I_{eq} = I_D(V_0) - g_m V_{GS0} - g_{ds} V_{DS0} - g_{mb} V_{BS0}$
 *
 * The KCL equations for the MOSFET (NMOS example):
 * Node D: $+I_D = ...$
 * Node G: $0 = ...$ (ideal gate, no current)
 * Node S: $-I_D = ...$
 * Node B: $0 = ...$ (ideal bulk, no current, unless bulk current is modeled)
 *
 * The current $I_D$ flows from Drain to Source.
 *
 * @param A The MNA matrix to which stamps are applied.
 * @param b The MNA right-hand side vector to which stamps are applied.
 * @param x_current_guess The current guess for node voltages and branch currents.
 * @param prev_solution The solution from the previous time step (not directly used here).
 * @param time The current simulation time (not used).
 * @param dt The time step size (not used).
 */
void MOSFET::getStamps(
    Eigen::MatrixXd& A, Eigen::VectorXd& b,
    const Eigen::VectorXd& x_current_guess, const Eigen::VectorXd& prev_solution,
    double time, double dt
) {
    // Get the global indices for the MOSFET's nodes.
    int idx_d = getNodeIndex(drain_node_name);
    int idx_g = getNodeIndex(gate_node_name);
    int idx_s = getNodeIndex(source_node_name);
    int idx_b = getNodeIndex(bulk_node_name);

    // Ensure all node indices are valid
    if (idx_d == -1 || idx_g == -1 || idx_s == -1 || idx_b == -1) {
        std::cerr << "Error: Invalid node index for MOSFET " << name << ". Check node names." << std::endl;
        return;
    }

    // Get current guess voltages from x_current_guess
    // Handle ground node (index 0) which has voltage 0
    double Vd_curr = (idx_d == 0) ? 0.0 : x_current_guess(idx_d - 1); // Adjust for 0-based indexing if node 0 is ground
    double Vg_curr = (idx_g == 0) ? 0.0 : x_current_guess(idx_g - 1);
    double Vs_curr = (idx_s == 0) ? 0.0 : x_current_guess(idx_s - 1);
    double Vb_curr = (idx_b == 0) ? 0.0 : x_current_guess(idx_b - 1);

    // Calculate terminal voltages
    double Vgs_curr = Vg_curr - Vs_curr;
    double Vds_curr = Vd_curr - Vs_curr;
    double Vbs_curr = Vb_curr - Vs_curr;
    double Vsb_curr = Vs_curr - Vb_curr; // Source-Bulk voltage for Vth calculation

    // Calculate threshold voltage with body effect
    double Vth_curr = calculateThresholdVoltage(Vsb_curr);

    // Calculate drain current and small-signal conductances at current operating point
    double Id_curr = calculateDrainCurrent(Vgs_curr, Vds_curr, Vth_curr);
    double gds_val = calculateGDS(Vgs_curr, Vds_curr, Vth_curr);
    double gm_val = calculateGM(Vgs_curr, Vds_curr, Vth_curr);
    double gmb_val = calculateGMB(Vgs_curr, Vds_curr, Vbs_curr, Vth_curr);

    // Adjust for PMOS: current direction is reversed, and some conductances might be too.
    // For PMOS, Id flows from Source to Drain.
    // The SPICE model equations are typically written for NMOS, and then applied to PMOS
    // by swapping D/S, using |Vgs|, |Vds|, |Vth|, and then inverting current.
    // Here, we calculate Id as if it's NMOS, and then reverse the stamp if PMOS.
    if (type == "PMOS") {
        // Id is calculated based on Vgs, Vds, Vth (which are negative for PMOS).
        // The current Id is positive if it flows from Source to Drain.
        // Our calculateDrainCurrent returns a value that is positive for NMOS in typical operation.
        // For PMOS, if Vgs < Vth (e.g., -2V < -1V), it's ON.
        // If Vds < Vgs - Vth (e.g., -3V < -2V - (-1V) = -1V), it's triode.
        // The current will be negative according to our NMOS-like calculation,
        // so we use its absolute value for the source-drain current, and then stamp accordingly.
        Id_curr = -Id_curr; // Id flows from Source to Drain for PMOS
    }


    // Apply stamps to the MNA matrix A and vector b
    // The MOSFET stamps are:
    // KCL at Drain: +Id_curr - gds*Vds_curr - gm*Vgs_curr - gmb*Vbs_curr
    // KCL at Source: -Id_curr + gds*Vds_curr + gm*Vgs_curr + gmb*Vbs_curr
    // The companion model current source I_eq = Id_curr - gds*Vds_curr - gm*Vgs_curr - gmb*Vbs_curr
    // This I_eq is added to the RHS of the KCL equation for the Drain node, and subtracted from Source.

    // KCL at Drain node (idx_d)
    // Add gds * Vd
    A(idx_d, idx_d) += gds_val;
    // Add gm * Vg
    A(idx_d, idx_g) += gm_val;
    // Add gmb * Vb
    A(idx_d, idx_b) += gmb_val;
    // Subtract (gds + gm + gmb) * Vs (since Vds = Vd-Vs, Vgs = Vg-Vs, Vbs = Vb-Vs)
    A(idx_d, idx_s) -= (gds_val + gm_val + gmb_val);

    // KCL at Source node (idx_s)
    // Subtract gds * Vd
    A(idx_s, idx_d) -= gds_val;
    // Subtract gm * Vg
    A(idx_s, idx_g) -= gm_val;
    // Subtract gmb * Vb
    A(idx_s, idx_b) -= gmb_val;
    // Add (gds + gm + gmb) * Vs
    A(idx_s, idx_s) += (gds_val + gm_val + gmb_val);

    // Calculate the equivalent current source for the RHS vector b
    // I_eq = Id_curr - gds*Vds_curr - gm*Vgs_curr - gmb*Vbs_curr
    double I_eq = Id_curr - gds_val * Vds_curr - gm_val * Vgs_curr - gmb_val * Vbs_curr;

    // Add I_eq to the RHS of the Drain KCL equation
    b(idx_d) -= I_eq; // Current flows OUT of Drain for NMOS (or IN for PMOS if Id_curr is negative)
    // Add -I_eq to the RHS of the Source KCL equation
    b(idx_s) += I_eq; // Current flows IN to Source for NMOS (or OUT for PMOS if Id_curr is negative)

    // Note: Gate and Bulk nodes are typically high impedance (no current contribution)
    // unless gate leakage or junction diodes are modeled, which are not part of basic Level 1.
}

