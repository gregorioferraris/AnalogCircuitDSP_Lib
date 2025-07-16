// File: SpeakerDriver.cpp
#include "SpeakerDriver.h"
#include <algorithm> // For std::clamp (C++17) or manual clamping

/**
 * @brief Implementation of the SpeakerDriver class constructor.
 * @param name Name of the component.
 * @param input_node Name of the input node.
 * @param output_node Name of the output node.
 * @param ground_node Name of the ground node.
 * @param sample_rate The audio sample rate.
 * @param saturation_threshold Threshold for non-linear excursion limiting.
 * @param low_cut_freq Cutoff frequency for the high-pass filter.
 * @param high_cut_freq Cutoff frequency for the low-pass filter.
 * @param resonance_freq Center frequency for the resonance peak.
 * @param resonance_q Q factor for the resonance peak.
 */
SpeakerDriver::SpeakerDriver(const std::string& name, const std::string& input_node, const std::string& output_node,
                             const std::string& ground_node, double sample_rate,
                             double saturation_threshold,
                             double low_cut_freq, double high_cut_freq,
                             double resonance_freq, double resonance_q)
    : Component(name, input_node, output_node, ground_node),
      sample_rate_(sample_rate),
      saturation_threshold_(saturation_threshold),
      // Initialize Biquad filters
      // High-pass filter for low-end rolloff
      low_cut_filter_(BiquadFilter::FilterType::HighPass, low_cut_freq, 0.707, sample_rate),
      // Low-pass filter for high-end rolloff
      high_cut_filter_(BiquadFilter::FilterType::LowPass, high_cut_freq, 0.707, sample_rate),
      // Peaking filter for resonance
      resonance_filter_(BiquadFilter::FilterType::Peaking, resonance_freq, resonance_q, sample_rate)
{
    // Ensure parameters are within reasonable bounds
    if (sample_rate_ <= 0) sample_rate_ = 44100.0;
    if (saturation_threshold_ <= 0) saturation_threshold_ = 0.01; // Avoid division by zero or infinite gain
    if (saturation_threshold_ > 10.0) saturation_threshold_ = 10.0; // Prevent extreme values

    // Update filter parameters to ensure they are set correctly after construction
    setLowCutFreq(low_cut_freq);
    setHighCutFreq(high_cut_freq);
    setResonanceFreq(resonance_freq);
    setResonanceQ(resonance_q);
}

/**
 * @brief Processes a single audio sample through the speaker driver model.
 *
 * The processing chain is as follows:
 * 1.  **Frequency Shaping:** The input sample passes through a series of Biquad filters
 * to simulate the speaker's frequency response:
 * -   `low_cut_filter_`: High-pass filter to remove sub-bass frequencies.
 * -   `high_cut_filter_`: Low-pass filter to remove harsh high frequencies.
 * -   `resonance_filter_`: Peaking filter to simulate the speaker's resonant frequency.
 * 2.  **Saturation/Excursion Limiting:** A soft saturation function (hyperbolic tangent)
 * is applied to simulate the physical excursion limits of the speaker cone at high volumes.
 * The `saturation_threshold_` controls the point at which this non-linearity becomes noticeable.
 *
 * @param input_sample The input sample to process.
 * @return The processed sample, simulating the output of a speaker.
 */
double SpeakerDriver::process(double input_sample) {
    // 1. Apply frequency shaping
    double filtered_sample = low_cut_filter_.process(input_sample);
    filtered_sample = high_cut_filter_.process(filtered_sample);
    filtered_sample = resonance_filter_.process(filtered_sample);

    // 2. Apply soft saturation for excursion limiting
    double output_sample = applySaturation(filtered_sample, saturation_threshold_);

    return output_sample;
}

/**
 * @brief Applies a soft saturation/clipping effect using tanh.
 *
 * This function scales the input sample by the inverse of the threshold,
 * applies the hyperbolic tangent function, and then scales it back.
 * This results in a smooth compression effect as the sample approaches
 * the threshold, simulating the non-linear behavior of a speaker at its limits.
 *
 * @param sample The input sample.
 * @param threshold The threshold at which saturation begins.
 * @return The saturated sample.
 */
double SpeakerDriver::applySaturation(double sample, double threshold) {
    // Avoid division by zero if threshold is extremely small
    if (threshold == 0.0) return 0.0;

    // Scale the input by the inverse of the threshold
    // Apply tanh function for soft saturation
    // Scale back to the original range
    return threshold * std::tanh(sample / threshold);
}

// Setter implementations
void SpeakerDriver::setSaturationThreshold(double threshold) {
    if (threshold > 0) saturation_threshold_ = threshold;
}

void SpeakerDriver::setLowCutFreq(double freq) {
    if (freq > 0 && freq < sample_rate_ / 2.0) {
        low_cut_filter_.setCutoffFrequency(freq);
        low_cut_filter_.calculateCoefficients(); // Recalculate coefficients after changing freq
    }
}

void SpeakerDriver::setHighCutFreq(double freq) {
    if (freq > 0 && freq < sample_rate_ / 2.0) {
        high_cut_filter_.setCutoffFrequency(freq);
        high_cut_filter_.calculateCoefficients(); // Recalculate coefficients
    }
}

void SpeakerDriver::setResonanceFreq(double freq) {
    if (freq > 0 && freq < sample_rate_ / 2.0) {
        resonance_filter_.setCutoffFrequency(freq);
        resonance_filter_.calculateCoefficients(); // Recalculate coefficients
    }
}

void SpeakerDriver::setResonanceQ(double q) {
    if (q > 0) {
        resonance_filter_.setQFactor(q);
        resonance_filter_.calculateCoefficients(); // Recalculate coefficients
    }
}
