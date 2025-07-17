// analog_stages/AnalogOutputStage.h
#ifndef ANALOG_OUTPUT_STAGE_H
#define ANALOG_OUTPUT_STAGE_H

#include <string>
#include <vector>
#include <memory> // For std::shared_ptr

// --- Forward Declarations for Components ---
// These are assumed to exist or will be defined elsewhere in your library.
// They are needed to use pointers to these types without including their .h files here.
class Component; // Base class for all components
class Resistor;
class Capacitor;

/**
 * @class AnalogOutputStage
 * @brief Represents an analog output stage as a subcircuit.
 *
 * This stage typically includes an output resistance and an optional low-pass filter
 * to model the output frequency response and prevent aliasing.
 * It acts as a container for simpler components.
 */
class AnalogOutputStage {
public:
    /**
     * @brief Constructor for the AnalogOutputStage class.
     * @param name Unique name of the output stage.
     * @param input_node Name of the input node of the stage.
     * @param output_node Name of the output node of the stage.
     * @param ground_node Name of the ground node.
     * @param output_resistance Value of the output resistance in Ohms.
     * @param filter_C Value of the low-pass filter capacitor in Farads (0 for no filter).
     */
    AnalogOutputStage(const std::string& name,
                      const std::string& input_node,
                      const std::string& output_node,
                      const std::string& ground_node,
                      double output_resistance,
                      double filter_C);

    // Method to access components (to pass them to the solver)
    const std::vector<std::shared_ptr<Component>>& getComponents() const { return components_; }

    // Method to get all internal and external nodes managed by this stage
    const std::vector<std::string>& getNodes() const { return nodes_; }

    // Method to get the main input/output nodes
    std::string getInputNode() const { return input_nodes_[0]; } // Assume a single input
    std::string getOutputNode() const { return output_nodes_[0]; } // Assume a single output

private:
    std::string name_;
    std::vector<std::string> input_nodes_;
    std::vector<std::string> output_nodes_;
    std::vector<std::string> nodes_; // All nodes, including internal ones
    std::vector<std::shared_ptr<Component>> components_; // Use shared_ptr

    // Internal nodes specific to this stage
    std::string filter_node_; // Intermediate node for the RC filter
    std::string ground_node_; // The ground reference node

    // Helper method to add a node (and ensure it's unique)
    void addNode(const std::string& node_name);

    // Helper method to add a component
    template<typename T, typename... Args>
    void addComponent(Args&&... args) {
        // Create an instance of the component and add it to the vector
        components_.push_back(std::make_shared<T>(std::forward<Args>(args)...));
    }
};

#endif // ANALOG_OUTPUT_STAGE_H
