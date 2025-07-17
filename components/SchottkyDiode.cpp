// File: SchottkyDiode.cpp
#include "SchottkyDiode.h"
#include <algorithm> // For std::max (though not strictly used in final version, good to keep)
#include <cmath>     // For std::atan

/**
 * @brief Implementation of the SchottkyDiode class constructor.
 * @param name Name of the component.
 * @param input_node Name of the input node.
 * @param output_node Name of the output node.
 * @param ground_node Name of the ground node.
 * @param forward_voltage_threshold Forward voltage (Vf).
 * @param softness Softness factor for the conduction curve.
 */
SchottkyDiode::SchottkyDiode(const std::string& name, const std::string& input_node, const std::string& output_node,
                             const std::string& ground_node, double forward_voltage_threshold, double softness)
    : Component(name, input_node, output_node, ground_node),
      forward_voltage_threshold_(forward_voltage_threshold),
      softness_(softness)
{
    // Ensure softness is positive to avoid division by zero or NaN.
    if (softness_ <= 0) {
        softness_ = 1e-6; // Set a small positive default value
    }
}

/**
 * @brief Processes a single audio sample through the Schottky diode model.
 *
 * The model simulates a half-wave rectification with a soft "knee".
 *
 * 1.  **Half-wave Rectification:** If the input sample is negative, the output is zero (the diode blocks).
 * 2.  **Conduction Threshold:** For positive samples, the `forward_voltage_threshold_` is subtracted.
 * This simulates the forward voltage drop required for conduction.
 * 3.  **Soft "Knee":** An `atan` (arctangent) function is applied to the signal
 * above the threshold. This creates a smoother conduction curve compared to hard clipping,
 * simulating the gradual turn-on behavior of a real diode. The `softness_`
 * parameter controls the slope and softness of this curve.
 *
 * @param input_sample The input sample to process.
 * @return The processed (rectified and shaped) sample.
 */
double SchottkyDiode::process(double input_sample) {
    double output_signal = 0.0;

    // A Schottky diode primarily conducts in one direction.
    // For a half-wave rectification model, we block negative signals.
    if (input_sample > 0) {
        // Calculate the signal voltage above the forward voltage threshold (Vf)
        double signal_above_threshold = input_sample - forward_voltage_threshold_;

        // If the signal is actually above the threshold, apply the conduction curve
        if (signal_above_threshold > 0) {
            // Use std::atan to create a soft "knee" in the conduction curve.
            // The 'softness_' controls the slope and softness of the knee.
            // We multiply by 'softness_' to keep the amplitude proportional
            // to the input for small signals above the thresho
