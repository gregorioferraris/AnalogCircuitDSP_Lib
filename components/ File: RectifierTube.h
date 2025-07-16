// File: RectifierTube.h
#ifndef RECTIFIER_TUBE_H
#define RECTIFIER_TUBE_H

#include "Component.h" // Assumes Component.h defines the base class Component
#include <string>
#include <cmath>     // For std::fabs, std::atan

/**
 * @class RectifierTube
 * @brief Models a simplified vacuum tube rectifier effect.
 *
 * This class simulates the non-linear behavior of a tube rectifier,
 * including full-wave rectification, a voltage drop characteristic,
 * and a soft "knee" in its conduction curve. It's an audio effect model,
 * not a full circuit simulation.
 */
class RectifierTube : public Component {
public:
    /**
     * @brief Constructor for the RectifierTube class.
     * @param name Name of the component.
     * @param input_node Name of the input node.
     * @param output_node Name of the output node.
     * @param ground_node Name of the ground node.
     * @param voltage_drop_factor Factor representing the internal voltage drop of the tube.
     * Higher values mean more voltage drop. (e.g., 0.05 for 5% drop)
     * @param softness_factor Factor controlling the softness of the rectification "knee".
     * A smaller value makes the knee sharper.
     */
    RectifierTube(const std::string& name, const std::string& input_node, const std::string& output_node,
                  const std::string& ground_node, double voltage_drop_factor = 0.05, double softness_factor = 0.1);

    /**
     * @brief Processes a single audio sample through the rectifier tube model.
     * @param input_sample The input sample to process (AC signal).
     * @return The processed (rectified and shaped) sample (pulsating DC).
     */
    double process(double input_sample) override;

    // Getter methods
    double getVoltageDropFactor() const { return voltage_drop_factor_; }
    double getSoftnessFactor() const { return softness_factor_; }

private:
    double voltage_drop_factor_; ///< Factor for simulating internal voltage drop.
    double softness_factor_;     ///< Softness factor for the rectification curve.
};

#endif // RECTIFIER_TUBE_H
