// File: SchottkyDiode.h
#ifndef SCHOTTKY_DIODE_H
#define SCHOTTKY_DIODE_H

#include "Component.h" // Assumes Component.h defines the base class Component
#include <string>
#include <cmath>     // For std::atan

/**
 * @class SchottkyDiode
 * @brief Models a simplified Schottky diode effect for clipping or rectification.
 *
 * This class simulates the non-linear behavior of a Schottky diode,
 * characterized by a low forward voltage drop and fast switching.
 * The model implements a half-wave rectification with a soft "knee"
 * around the conduction threshold.
 */
class SchottkyDiode : public Component {
public:
    /**
     * @brief Constructor for the SchottkyDiode class.
     * @param name Name of the component.
     * @param input_node Name of the input node.
     * @param output_node Name of the output node.
     * @param ground_node Name of the ground node.
     * @param forward_voltage_threshold Forward voltage (Vf) threshold of the diode.
     * @param softness Factor controlling the softness of the conduction "knee".
     * A smaller value makes the knee sharper.
     */
    SchottkyDiode(const std::string& name, const std::string& input_node, const std::string& output_node,
                  const std::string& ground_node, double forward_voltage_threshold = 0.3, double softness = 0.05);

    /**
     * @brief Processes a single audio sample through the Schottky diode model.
     * @param input_sample The input sample to process.
     * @return The processed (rectified and shaped) sample.
     */
    double process(double input_sample) override;

    // Getter methods
    double getForwardVoltageThreshold() const { return forward_voltage_threshold_; }
    double getSoftness() const { return softness_; }

private:
    double forward_voltage_threshold_; ///< Forward voltage threshold of the diode.
    double softness_;                  ///< Softness factor for the conduction curve.
};

#endif // SCHOTTKY_DIODE_H
