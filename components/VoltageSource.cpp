// File: VoltageSource.cpp
#include "VoltageSource.h"
#include <cmath> // For M_PI (or define PI if not available)

// Define PI if M_PI is not available in cmath on some compilers
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/**
 * @brief Implementation of the VoltageSource class constructor.
 * @param name Name of the component.
 * @param output_node Name of the output node.
 * @param ground_node Name of the ground node.
 * @param type The type of waveform to generate.
 * @param initial_voltage Initial DC voltage or amplitude for AC waveforms.
 * @param frequency Frequency for AC waveforms (Hz).
 * @param phase_offset Phase offset for AC waveforms (radians).
 * @param sample_rate The audio sample rate in Hz.
 */
VoltageSource::VoltageSource(const std::string& name, const std::string& output_node,
                             const std::string& ground_node, WaveformType type,
                             double initial_voltage, double frequency,
                             double phase_offset, double sample_rate)
    : Component(name, "", output_node, ground_node), // VoltageSource typically has no input node itself, it generates output
      type_(type),
      voltage_(initial_voltage),
      frequency_(frequency),
      phase_offset_(phase_offset),
      sample_rate_(sample_rate),
      current_phase_(0.0)
{
    // Ensure sample rate is valid for AC waveforms
    if (sample_rate_ <= 0) sample_rate_ = 44100.0;
    // Ensure frequency is non-negative
    if (frequency_ < 0) frequency_ = 0.0;
}

/**
 * @brief Processes a single sample for the voltage source.
 *
 * This method generates the appropriate voltage based on the `WaveformType` configured.
 * - `DC`: Returns the constant `voltage_`.
 * - `SINE`: Generates a sine wave based on `voltage_` (amplitude), `frequency_`,
 * `phase_offset_`, and `sample_rate_`. The `current_phase_` is updated for
 * the next sample.
 * - `SQUARE`: Generates a square wave. The output is `voltage_` for the first
 * half of the period and `-voltage_` for the second half.
 * - `EXTERNAL`: Simply returns the `input_sample` provided, acting as a pass-through.
 *
 * @param input_sample The input sample (only used for EXTERNAL type).
 * @return The generated or passed-through voltage.
 */
double VoltageSource::process(double input_sample) {
    double output_voltage = 0.0;

    switch (type_) {
        case DC:
            output_voltage = voltage_;
            break;
        case SINE: {
            // Calculate instantaneous phase in radians
            double phase_radians = current_phase_ + phase_offset_;
            output_voltage = voltage_ * std::sin(phase_radians);

            // Update phase for the next sample
            current_phase_ += (2.0 * M_PI * frequency_) / sample_rate_;
            // Keep phase within 0 to 2*PI to prevent large numbers and precision issues
            current_phase_ = std::fmod(current_phase_, 2.0 * M_PI);
            break;
        }
        case SQUARE: {
            // Calculate instantaneous phase in radians for square wave logic
            double phase_radians = current_phase_ + phase_offset_;
            // Normalize phase to 0-1 range for easier square wave logic
            double normalized_phase = std::fmod(phase_radians / (2.0 * M_PI), 1.0);
            if (normalized_phase < 0) normalized_phase += 1.0; // Ensure positive

            if (normalized_phase < 0.5) { // First half of the period
                output_voltage = voltage_;
            } else { // Second half of the period
                output_voltage = -voltage_;
            }

            // Update phase for the next sample
            current_phase_ += (2.0 * M_PI * frequency_) / sample_rate_;
            current_phase_ = std::fmod(current_phase_, 2.0 * M_PI);
            break;
        }
        case EXTERNAL:
            output_voltage = input_sample;
            break;
        default:
            // Should not happen, but as a fallback
            output_voltage = 0.0;
            break;
    }

    return output_voltage;
}

// Setter implementations
void VoltageSource::setWaveformType(WaveformType type) {
    type_ = type;
    // Reset phase if changing to an AC type to ensure a consistent start
    if (type == SINE || type == SQUARE) {
        current_phase_ = 0.0;
    }
}

void VoltageSource::setVoltage(double voltage) {
    voltage_ = voltage;
}

void VoltageSource::setFrequency(double freq) {
    if (freq >= 0) { // Frequency cannot be negative
        frequency_ = freq;
        // Reset phase if frequency changes to avoid jumps in waveform
        current_phase_ = 0.0;
    }
}

void VoltageSource::setPhaseOffset(double offset) {
    phase_offset_ = offset;
}

void VoltageSource::setSampleRate(double rate) {
    if (rate > 0) { // Sample rate must be positive
        sample_rate_ = rate;
        // Reset phase if sample rate changes to avoid jumps in waveform
        current_phase_ = 0.0;
    }
}
