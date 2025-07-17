// File: RectifierTube.cpp
#include "RectifierTube.h"
#include <algorithm> // For std::max (though not strictly used in final version, good to keep)
#include <cmath> // For std::fabs, std::atan

/**
 * @brief Implementation of the RectifierTube class constructor.
 * @param name Name of the component.
 * @param input_node Name of the input node.
 * @param output_node Name of the output node.
 * @param ground_node Name of the ground node.
 * @param voltage_drop_factor Factor for internal voltage drop.
 * @param softness_factor Softness factor for the rectification curve.
 */
RectifierTube::RectifierTube(const std::string& name, const std::string& input_node, const std::string& output_node,
                             const std::string& ground_node, double voltage_drop_factor, double softness_factor)
    : Component(name, input_node, output_node, ground_node),
      voltage_drop_factor_(voltage_drop_factor),
      softness_factor_(softness_factor)
{
    // Ensure factors are within reasonable bounds to prevent unexpected behavior.
    if (voltage_drop_factor_ < 0) {
        voltage_drop_factor_ = 0.0;
    } else if (voltage_drop_factor_ > 1.0) { // Max 100% drop
        voltage_drop_factor_ = 1.0;
    }

    if (softness_factor_ <= 0) {
        softness_factor_ = 1e-6; // Set a small positive default value
    }
}

/**
 * @brief Processes a single audio sample through the rectifier tube model.
 *
 * The model simulates a full-wave rectification with a non-linear voltage drop
 * and a soft "knee" characteristic of vacuum tube rectifiers.
 *
 * 1.  **Full-wave Rectification:** The absolute value of the input sample is taken.
 * This converts both positive and negative parts of the AC signal into positive values.
 * 2.  **Simulated Voltage Drop:** A voltage drop is applied, proportional to the absolute
 * input signal, controlled by `voltage_drop_factor_`. This simulates the internal
 * resistance and voltage "sag" of a tube rectifier.
 * 3.  **Soft "Knee" Conduction:** An `atan` (arctangent) function is applied to the
 * signal after the voltage drop. This creates a smoother, more gradual turn-on
 * characteristic (a "soft knee") compared to ideal diode rectification,
 * which is typical of vacuum tubes. The `softness_factor_` controls
 * the slope and the degree of this softness.
 *
 * @param input_sample The input sample to process.
 * @return The processed (rectified and shaped) sample.
 */
double RectifierTube::process(double input_sample) {
    // 1. Full-wave rectification: Take the absolute value of the input.
    double rectified_signal = std::fabs(input_sample);

    // 2. Simulate internal voltage drop.
    // This simple model subtracts a percentage of the rectified signal,
    // simulating a non-ideal voltage drop under load.
    double signal_after_drop = rectified_signal * (1.0 - voltage_drop_factor_);

    // Ensure signal doesn't go negative after drop, though atan handles this well.
    if (signal_after_drop < 0) {
        signal_after_drop = 0;
    }

    // 3. Apply soft "knee" using atan.
    // This simulates the non-linear, gradual turn-on characteristic of a tube.
    // The softness_factor_ scales the input to atan and the output, controlling the curve.
    double final_output = softness_factor_ * std::atan(signal_after_drop / softness_factor_);

    return final_output;
}

/**
 * @brief Updates the internal state of the RectifierTube.
 *
 * This component is a static audio effect, so it does not have
 * internal state that evolves over time. This method is empty.
 *
 * @param v_curr The current voltage (not used).
 * @param i_curr The current (not used).
 */
void RectifierTube::updateState(double v_curr, double i_curr) {
    // No state to update for this static audio effect.
}
