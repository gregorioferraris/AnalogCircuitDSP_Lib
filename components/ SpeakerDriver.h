// File: SpeakerDriver.h
#ifndef SPEAKER_DRIVER_H
#define SPEAKER_DRIVER_H

#include "Component.h" // Assumes Component.h defines the base class Component
#include "BiquadFilter.h" // Assumes BiquadFilter.h defines a biquad filter class
#include <string>
#include <cmath>     // For std::tanh (hyperbolic tangent for saturation)

/**
 * @class SpeakerDriver
 * @brief Models a simplified speaker driver effect.
 *
 * This class simulates some non-linear and frequency-dependent characteristics
 * of a loudspeaker, including excursion limiting (saturation) and a basic
 * frequency response shaping. It's an audio effect model, not a full
 * electro-mechanical simulation.
 */
class SpeakerDriver : public Component {
public:
    /**
     * @brief Constructor for the SpeakerDriver class.
     * @param name Name of the component.
     * @param input_node Name of the input node (e.g., "amplifier_output").
     * @param output_node Name of the output node (e.g., "mic_input").
     * @param ground_node Name of the ground node.
     * @param sample_rate The audio sample rate in Hz.
     * @param saturation_threshold Threshold for non-linear excursion limiting (e.g., 0.8 for 80% of max).
     * @param low_cut_freq Cutoff frequency for the high-pass filter (Hz).
     * @param high_cut_freq Cutoff frequency for the low-pass filter (Hz).
     * @param resonance_freq Center frequency for the resonance peak (Hz).
     * @param resonance_q Q factor for the resonance peak.
     */
    SpeakerDriver(const std::string& name, const std::string& input_node, const std::string& output_node,
                  const std::string& ground_node, double sample_rate,
                  double saturation_threshold = 1.0,
                  double low_cut_freq = 80.0, double high_cut_freq = 5000.0,
                  double resonance_freq = 200.0, double resonance_q = 0.707);

    /**
     * @brief Processes a single audio sample through the speaker driver model.
     * @param input_sample The input sample to process (e.g., amplifier output).
     * @return The processed sample (simulated speaker output).
     */
    double process(double input_sample) override;

    // Setter methods for parameters (optional, but good for real-time adjustments)
    void setSaturationThreshold(double threshold);
    void setLowCutFreq(double freq);
    void setHighCutFreq(double freq);
    void setResonanceFreq(double freq);
    void setResonanceQ(double q);

    // Getter methods
    double getSaturationThreshold() const { return saturation_threshold_; }
    double getLowCutFreq() const { return low_cut_filter_.getCutoffFrequency(); }
    double getHighCutFreq() const { return high_cut_filter_.getCutoffFrequency(); }
    double getResonanceFreq() const { return resonance_filter_.getCutoffFrequency(); }
    double getResonanceQ() const { return resonance_filter_.getQFactor(); }

private:
    double sample_rate_;            ///< Audio sample rate.
    double saturation_threshold_;   ///< Threshold for non-linear saturation.

    BiquadFilter low_cut_filter_;   ///< High-pass filter for low-end rolloff.
    BiquadFilter high_cut_filter_;  ///< Low-pass filter for high-end rolloff.
    BiquadFilter resonance_filter_; ///< Peaking/Band-pass filter for resonance.

    /**
     * @brief Applies a soft saturation/clipping effect using tanh.
     * @param sample The input sample.
     * @param threshold The threshold at which saturation begins.
     * @return The saturated sample.
     */
    double applySaturation(double sample, double threshold);
};

#endif // SPEAKER_DRIVER_H
