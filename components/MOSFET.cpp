// components/PNP_BJT.cpp
#include "PNP_BJT.h"
#include <iostream> // For error messages
#include <cmath>    // For std::exp, std::fabs
#include <algorithm> // For std::max, std::min

// Define a small epsilon for numerical stability
const double BJT_EPSILON = 1e-9; // Small voltage difference for numerical stability

/**
 * @brief Constructor for the PNP_BJT component.
 *
 * Initializes a PNP BJT with its name, three connected nodes (Base, Collector, Emitter),
 * and key parameters for the Ebers-Moll model.
 *
 * @param name The unique name of the BJT.
 * @param node_base The name of the Base node.
 * @param node_collector The name of the Collector node.
 * @param node_emitter The name of the Emitter node.
 * @param Is The saturation current (Amps). Must be > 0.
 * @param beta_F The forward common-emitter current gain (Beta_F). Must be > 0.
 * @param Vt The thermal voltage (Volts). Default is 0.02585V (approx at 300K). Must be > 0.
 * @param ideality_factor The ideality factor (n). Must be > 0.
 * @param tolerance_percent Optional tolerance percentage.
 */
PNP_BJT::PNP_BJT(const std::string& name,
                 const std::string& node_base, const std::string& node_collector, const std::string& node_emitter,
                 double Is, double beta_F, double Vt, double ideality_factor, double tolerance_percent)
    // The base Component constructor takes node1 and node2.
    // For a 3-terminal device, we can map node1 to Base and node2 to Emitter for the base class,
    // and handle node_collector separately. The base Component's node1 and node2
    // are not strictly used for the BJT's internal connections, but for its
    // overall presence in the circuit's node list.
    : Component(name, node_base, node_emitter, tolerance_percent), // node1=base, node2=emitter
      node_base_name(node_base), node_collector_name(node_collector), node_emitter_name(node_emitter),
      Is_(Is), beta_F_(beta_F), Vt_(Vt), ideality_factor_(ideality_factor)
{
    // A BJT introduces two auxiliary variables for its internal diode currents (Ibe, Ibc)
    // or one for the collector current if using a simplified model.
    // For a full Ebers-Moll companion model, we often use controlled sources,
    // which don't necessarily require auxiliary variables if stamped directly.
    // However, if modeling internal diodes as separate branches, they might.
    // Let's assume for now that the stamping will be direct using conductances and current sources.
    setNumAuxiliaryVariables(0); // No explicit auxiliary variables for this model (currents are derived)

    // Add the collector node to the list of nodes managed by the component.
    // This is crucial for getNodeIndex to work correctly for the collector node.
    addNode(node_collector);

    // Ensure parameters are reasonable
    if (Is_ <= 0) Is_ = 1e-15; // Small positive value
    if (beta_F_ <= 0) beta_F_ = 100.0; // Default beta
    if (Vt_ <= 0) Vt_ = 0.02585; // Thermal voltage at 300K
    if (ideality_factor_ <= 0) ideality_factor_ = 1.0;

    std::cout << "PNP_BJT " << name << " initialized. B:" << node_base << " C:" << node_collector << " E:" << node_emitter << "." << std::endl;
}

/**
 * @brief Placeholder for processing a single audio sample (not used in MNA stamping).
 *
 * This method is inherited from Component but is not directly used when the BJT
 * is integrated via MNA stamps. Its purpose would be for a simpler "effect" model.
 *
 * @param input_sample The input sample.
 * @return A dummy value (0.0) as this model is for MNA.
 */
double PNP_BJT::process(double input_sample) {
    // This method is not used for MNA stamping.
    // The actual current calculation happens in getStamps based on node voltages.
    return 0.0;
}

/**
 * @brief Calculates the diode current for a PN junction.
 * @param V The voltage across the diode.
 * @param Is The saturation current.
 * @param Vt The thermal voltage.
 * @param n The ideality factor.
 * @return The diode current.
 */
double PNP_BJT::calculateDiodeCurrent(double V, double Is, double Vt, double n) const {
    double exp_arg = V / (n * Vt);

    // Clamp the exponential argument to prevent numerical overflow/underflow.
    // For very large positive V, exp_arg can be huge. For very large negative V, it can be very small.
    const double MAX_EXP_ARG = 30.0; // Corresponds to V ~ 0.77V for n=1, Vt=0.02585V
    const double MIN_EXP_ARG = -30.0; // Corresponds to V ~ -0.77V for n=1, Vt=0.02585V

    if (exp_arg > MAX_EXP_ARG) {
        exp_arg = MAX_EXP_ARG;
    } else if (exp_arg < MIN_EXP_ARG) {
        exp_arg = MIN_EXP_ARG;
    }

    return Is * (std::exp(exp_arg) - 1.0);
}

/**
 * @brief Calculates the small-signal conductance (transconductance) of a diode.
 * g = dI/dV = (Is / (n*Vt)) * exp(V / (n*Vt))
 * @param V The voltage across the diode.
 * @param Is The saturation current.
 * @param Vt The thermal voltage.
 * @param n The ideality factor.
 * @return The diode conductance.
 */
double PNP_BJT::calculateDiodeConductance(double V, double Is, double Vt, double n) const {
    double exp_arg = V / (n * Vt);

    // Clamp the exponential argument for derivative as well
    const double MAX_EXP_ARG = 30.0;
    const double MIN_EXP_ARG = -30.0;

    if (exp_arg > MAX_EXP_ARG) {
        exp_arg = MAX_EXP_ARG;
    } else if (exp_arg < MIN_EXP_ARG) {
        exp_arg = MIN_EXP_ARG;
    }

    return (Is / (n * Vt)) * std::exp(exp_arg);
}


/**
 * @brief Applies the stamps of the PNP BJT companion model to the MNA matrix (A) and vector (b).
 *
 * This method implements a simplified Ebers-Moll model for the PNP BJT,
 * linearizing it for MNA using conductances and equivalent current sources.
 *
 * For a PNP BJT:
 * - Base-Emitter (BE) junction acts like a diode. Vbe = Vb - Ve.
 * - Base-Collector (BC) junction acts like a diode. Vbc = Vb - Vc.
 *
 * Currents (simplified Ebers-Moll):
 * I_C = -alpha_F * I_F + I_R (Collector current, positive into collector)
 * I_E = I_F - alpha_R * I_R (Emitter current, positive into emitter)
 * I_B = I_F * (1 - alpha_F) + I_R * (1 - alpha_R) (Base current, positive into base)
 *
 * Where:
 * I_F = Is * (exp(Vbe / Vt) - 1) (Forward current due to BE diode)
 * I_R = Is * (exp(Vbc / Vt) - 1) (Reverse current due to BC diode)
 * alpha_F = Beta_F / (1 + Beta_F) (Forward common-base current gain)
 * alpha_R = Beta_R / (1 + Beta_R) (Reverse common-base current gain, often assumed small or zero for simplified model)
 *
 * For a simplified model, often only the forward active region is considered:
 * I_C = -beta_F * I_B_diode (current flowing out of collector)
 * I_B = I_B_diode
 * where I_B_diode is the current through the BE junction.
 *
 * Let's use the standard companion model approach for the two diodes (BE and BC)
 * and then add the controlled source.
 *
 * Diode model: I = G_eq * V + I_eq
 * G_eq = dI/dV
 * I_eq = I_actual - G_eq * V_actual
 *
 * For PNP:
 * Vbe_curr = Vb_curr - Ve_curr
 * Vbc_curr = Vb_curr - Vc_curr
 *
 * Calculate currents and conductances for BE and BC diodes:
 * I_BE = Is_ * (exp(Veb_curr / (ideality_factor_ * Vt_)) - 1) where Veb_curr = Ve_curr - Vb_curr
 * G_BE = (Is_ / (ideality_factor_ * Vt_)) * exp(Veb_curr / (ideality_factor_ * Vt_))
 *
 * I_BC = Is_ * (exp(Vcb_curr / (ideality_factor_ * Vt_)) - 1) where Vcb_curr = Vc_curr - Vb_curr
 * G_BC = (Is_ / (ideality_factor_ * Vt_)) * exp(Vcb_curr / (ideality_factor_ * Vt_))
 *
 * Equivalent current sources:
 * I_BE_eq = I_BE - G_BE * Veb_curr
 * I_BC_eq = I_BC - G_BC * Vcb_curr
 *
 * Transconductance gm for PNP (current from Collector to Emitter controlled by Veb):
 * gm = d(Ic)/d(Veb) = beta_F * d(Ibe)/d(Veb) = beta_F * G_BE
 *
 * Now, apply stamps:
 * KCL at Base (idx_b):
 * - From BE diode (current I_BE_diode flows E->B, so it flows INTO Base from Emitter)
 * - From BC diode (current I_BC_diode flows C->B, so it flows INTO Base from Collector)
 *
 * KCL at Collector (idx_c):
 * - From BC diode (current I_BC_diode flows C->B, so it flows OUT of Collector to Base)
 * - From controlled current source (I_C_controlled_actual flows E->C, so it flows INTO Collector from Emitter)
 *
 * KCL at Emitter (idx_e):
 * - From BE diode (current I_BE_diode flows E->B, so it flows OUT of Emitter to Base)
 * - From controlled current source (I_C_controlled_actual flows E->C, so it flows OUT of Emitter to Collector)
 *
 * @param A The MNA matrix to which stamps are applied.
 * @param b The MNA right-hand side vector to which stamps are applied.
 * @param x_current_guess The current guess for node voltages and branch currents.
 * @param prev_solution The solution from the previous time step (not directly used here for static non-linearity).
 * @param time The current simulation time (not used).
 * @param dt The time step size (not used for static non-linearity).
 */
void PNP_BJT::getStamps(
    Eigen::MatrixXd& A, Eigen::VectorXd& b,
    const Eigen::VectorXd& x_current_guess, const Eigen::VectorXd& prev_solution,
    double time, double dt
) {
    // Get the global (1-based) indices for the BJT's nodes.
    int idx_b = getNodeIndex(node_base_name);
    int idx_c = getNodeIndex(node_collector_name);
    int idx_e = getNodeIndex(node_emitter_name);

    // Ensure all node indices are valid
    if (idx_b == -1 || idx_c == -1 || idx_e == -1) {
        std::cerr << "Error: Invalid node index for PNP_BJT " << name << ". Check node names." << std::endl;
        return;
    }

    // Get current guess voltages from x_current_guess
    // Adjust for 0-based indexing if node 0 is ground
    double Vb_curr = (idx_b == 0) ? 0.0 : x_current_guess(idx_b - 1);
    double Vc_curr = (idx_c == 0) ? 0.0 : x_current_guess(idx_c - 1);
    double Ve_curr = (idx_e == 0) ? 0.0 : x_current_guess(idx_e - 1);

    // Calculate terminal voltages
    double Veb_curr = Ve_curr - Vb_curr; // Emitter-Base voltage (positive for forward bias of BE diode)
    double Vcb_curr = Vc_curr - Vb_curr; // Collector-Base voltage (positive for forward bias of BC diode)

    // --- BE Diode (Emitter-Base) ---
    // Current I_BE_diode flows from Emitter to Base
    double I_BE_diode = calculateDiodeCurrent(Veb_curr, Is_, Vt_, ideality_factor_);
    double G_BE_diode = calculateDiodeConductance(Veb_curr, Is_, Vt_, ideality_factor_);

    // Equivalent current source for BE diode: I_eq = I_actual - G_eq * V_actual
    double I_BE_eq = I_BE_diode - G_BE_diode * Veb_curr;

    // --- BC Diode (Collector-Base) ---
    // Current I_BC_diode flows from Collector to Base
    double I_BC_diode = calculateDiodeCurrent(Vcb_curr, Is_, Vt_, ideality_factor_);
    double G_BC_diode = calculateDiodeConductance(Vcb_curr, Is_, Vt_, ideality_factor_);

    // Equivalent current source for BC diode
    double I_BC_eq = I_BC_diode - G_BC_diode * Vcb_curr;

    // --- Controlled Current Source (Collector Current) ---
    // In forward active region, Ic = beta_F * I_BE_diode (current from Emitter to Collector)
    // This current is controlled by Veb.
    // I_C_controlled_actual = beta_F_ * I_BE_diode
    // Transconductance gm = d(I_C_controlled_actual)/d(Veb) = beta_F_ * G_BE_diode

    double I_C_controlled_actual = beta_F_ * I_BE_diode; // Current from Emitter to Collector
    double gm_val = beta_F_ * G_BE_diode; // Transconductance from Veb to Ic

    // Equivalent current source for the controlled collector current
    double I_C_controlled_eq = I_C_controlled_actual - gm_val * Veb_curr;

    // --- Apply MNA stamps ---

    // KCL at Base node (idx_b)
    // From BE diode (current I_BE_diode flows E->B, so it flows INTO Base from Emitter)
    // I_B = I_BE_diode + I_BC_diode
    // Stamp for I_BE_diode:
    A(idx_b, idx_b) += G_BE_diode; // Vb coefficient
    A(idx_b, idx_e) -= G_BE_diode; // Ve coefficient
    b(idx_b) -= I_BE_eq;           // Equivalent current source (current flows INTO node)

    // Stamp for I_BC_diode:
    A(idx_b, idx_b) += G_BC_diode; // Vb coefficient
    A(idx_b, idx_c) -= G_BC_diode; // Vc coefficient
    b(idx_b) -= I_BC_eq;           // Equivalent current source (current flows INTO node)


    // KCL at Collector node (idx_c)
    // From BC diode (current I_BC_diode flows C->B, so it flows OUT of Collector to Base)
    // Stamp for I_BC_diode:
    A(idx_c, idx_c) += G_BC_diode; // Vc coefficient
    A(idx_c, idx_b) -= G_BC_diode; // Vb coefficient
    b(idx_c) += I_BC_eq;           // Equivalent current source (current flows OUT of node)

    // From controlled current source (I_C_controlled_actual flows E->C, so it flows INTO Collector from Emitter)
    // Stamp for I_C_controlled_actual:
    A(idx_c, idx_e) += gm_val; // Ve coefficient (for Veb = Ve - Vb, so gm * Ve)
    A(idx_c, idx_b) -= gm_val; // Vb coefficient (for Veb = Ve - Vb, so -gm * Vb)
    b(idx_c) -= I_C_controlled_eq; // Equivalent current source (current flows INTO node)


    // KCL at Emitter node (idx_e)
    // From BE diode (current I_BE_diode flows E->B, so it flows OUT of Emitter to Base)
    // Stamp for I_BE_diode:
    A(idx_e, idx_e) += G_BE_diode; // Ve coefficient
    A(idx_e, idx_b) -= G_BE_diode; // Vb coefficient
    b(idx_e) += I_BE_eq;           // Equivalent current source (current flows OUT of node)

    // From controlled current source (I_C_controlled_actual flows E->C, so it flows OUT of Emitter to Collector)
    // Stamp for I_C_controlled_actual:
    A(idx_e, idx_b) += gm_val; // Vb coefficient (for Veb = Ve - Vb, so gm * (-Vb))
    A(idx_e, idx_e) -= gm_val; // Ve coefficient (for Veb = Ve - Vb, so -gm * Ve)
    b(idx_e) += I_C_controlled_eq; // Equivalent current source (current flows OUT of node)
}

/**
 * @brief Updates the internal state of the PNP BJT.
 *
 * For a static non-linear model, there is no internal state to update
 * based on previous time steps. This method is empty.
 *
 * @param v_curr The current voltage across the component (not used).
 * @param i_curr The current flowing through the component (not used).
 */
void PNP_BJT::updateState(double v_curr, double i_curr) {
    // No state to update for this static non-linear model.
}
