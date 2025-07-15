// circuit_solver/Circuit.cpp
#include "Circuit.h"

Circuit::Circuit() : nextNodeId(1) {
    nodes["0"] = 0; // Node '0' is always ground and has ID 0
}

int Circuit::addNode(const std::string& node_name) {
    if (nodes.find(node_name) == nodes.end()) {
        nodes[node_name] = nextNodeId;
        nextNodeId++;
    }
    return nodes[node_name];
}

void Circuit::addComponent(std::shared_ptr<Component> component) {
    // Assign a unique ID to the component
    component->setComponentId(static_cast<int>(components.size()));

    // Assign numeric IDs to the component's nodes
    std::vector<int> node_ids_for_component;
    for (const std::string& node_name : component->getNodeNamesStr()) {
        node_ids_for_component.push_back(addNode(node_name));
    }
    component->setNodeIds(node_ids_for_component);

    // Add the component to the general list
    components.push_back(component);

    // Categorize components that add auxiliary variables
    // dynamic_pointer_cast is safe for downcasting with shared_ptr
    if (std::shared_ptr<VoltageSource> vs = std::dynamic_pointer_cast<VoltageSource>(component)) {
        voltageSources.push_back(vs);
    }
    if (std::shared_ptr<Splitter> splitter = std::dynamic_pointer_cast<Splitter>(component)) {
        splitters.push_back(splitter);
    }

    std::cout << "Added component: " << component->getName() << " (" << typeid(*component).name() << ")" << std::endl;
}

// The implementation of connectBlocks() would be here if we were using Subcircuit in C++
// For now, the management of connections between MNA blocks and functional blocks
// will be handled at the main.cpp level, as in your last Python example.
