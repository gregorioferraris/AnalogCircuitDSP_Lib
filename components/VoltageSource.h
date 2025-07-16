// File: VoltageSource.h
#ifndef VOLTAGE_SOURCE_H
#define VOLTAGE_SOURCE_H

#include "Component.h" // Assumes Component.h defines the base class Component
#include <string>
#include <cmath>     // For std::sin, std::fmod

/**
 * @class VoltageSource
 * @brief Models a simplified voltage source component.
 *
 * This class can generate various types of voltage signals, including
 * DC (constant), sine wave, square wave, or pass through an external value.
 * It's a fundamental component for providing input signals to a circuit simulation.
 */
class VoltageSource : public Component {
public:
    /**
     * @brief Enum for different types of voltage waveforms.
     */
    enum WaveformType {
        DC,         ///< Direct Current (constant voltage)
        SINE,       ///< Sine wave
        SQUARE,     ///< Square wave
        EXTERNAL    ///< External input (passes through the provided sample)
    };

    /**
     * @brief Constructor for the VoltageSource class.
     * @param name Name of the component.
     * @param output_node Name of the output node (where the voltage is produced).
     * @param ground_node Name of the ground node.
     * @param type The type of waveform to generate.
     * @param initial_voltage Initial DC voltage or amplitude for AC waveforms.
     * @param frequency Frequency for AC waveforms (Hz).
     * @param phase_offset Phase offset for AC waveforms (radians).
     * @param sample_rate The audio sample rate in Hz (needed for AC waveforms).
     */
    VoltageSource(const std::string& name, const std::string& output_node,
                  const std::string& ground_node, WaveformType type = DC,
                  double initial_voltage = 0.0, double frequency = 1.0,
                  double phase_offset = 0.0, double sample_rate = 44100.0);

    /**
     * @brief Processes a single sample for the voltage source.
     *
     * For DC, SINE, and SQUARE types, it generates the voltage.
     * For EXTERNAL type, it simply returns the input_sample.
     *
     * @param input_sample The input sample (only used for EXTERNAL type).
     * @return The generated or passed-through voltage.
     */
    double process(double input_sample) override;

    // Setter methods for parameters
    void setWaveformType(WaveformType type);
    void setVoltage(double voltage); // For DC or amplitude for AC
    void setFrequency(double freq);   // For AC waveforms
    void setPhaseOffset(double offset); // For AC waveforms
    void setSampleRate(double rate);    // For AC waveforms

    // Getter methods
    WaveformType getWaveformType() const { return type_; }
    double getVoltage() const { return voltage_; }
    double getFrequency() const { return frequency_; }
    double getPhaseOffset() const { return phase_offset_; }
    double getSampleRate() const { return sample_rate_; }

private:
    WaveformType type_;         ///< Type of waveform to generate.
    double voltage_;            ///< DC voltage or amplitude for AC.
    double frequency_;          ///< Frequency for AC waveforms (Hz).
    double phase_offset_;       ///< Phase offset for AC waveforms (radians).
    double sample_rate_;        ///< Audio sample rate (Hz).
    double current_phase_;      ///< Current phase for AC waveform generation.
};

#endif // VOLTAGE_SOURCE_H
